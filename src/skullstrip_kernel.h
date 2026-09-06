// Sample the working volume by NEAREST NEIGHBOUR, never by interpolation. AFNI's ray walk
// does the same, and it matters: interpolating smooths the ray minimum away -- measured on
// T1w node 1000 at it=0, trilinear gave Imin 64.6 where the reference reports 59.
// The volume is indexed in voxel units, so origin 0 and delta 1.
static float SS_FN(ss_sample_nn)(const unsigned char *vol, int nx, int ny, int nz,
		float x, float y, float z, int *out) {
	// IN-RANGE FAST PATH, in BOTH specialisations: the origin-0 delta-1 case of
	// ss_world_to_index with all three axes tested by ONE unsigned compare each -- a negative
	// index wraps to a huge unsigned, so `< n` rejects both ends at once and no clamp is
	// needed. The ray walk STOPS on the first out-of-range step, so at most one step per ray
	// ever needs the clamp; making every step pay six signed compare/select pairs to serve it
	// was the wrong trade (measured: 11% of total op runtime). ss_world_to_index remains the
	// definition and still runs whenever the answer is not simply `i` -- do not duplicate its
	// clamp here, and if its conversion convention ever changes, this test must follow it.
	//
	// Where the two specialisations differ is only WHERE the cold half lives. Inline, clang
	// if-converts it and computes the clamp on every step anyway; out of line it does not.
	// That is worth a further ~8%, but it also moves the result, so faithful keeps it inline.
	int i = (int)(x + 0.49f), j = (int)(y + 0.49f), k = (int)(z + 0.49f);

	if ((unsigned)i < (unsigned)nx && (unsigned)j < (unsigned)ny &&
			(unsigned)k < (unsigned)nz) {
		*out = 0;
	} else {
#ifdef SS_FAST
		ss_sample_clamp(nx, ny, nz, x, y, z, &i, &j, &k, out);
#else
		int oi, oj, ok;
		i = ss_world_to_index(x, 0.0f, 1.0f, nx, &oi);
		j = ss_world_to_index(y, 0.0f, 1.0f, ny, &oj);
		k = ss_world_to_index(z, 0.0f, 1.0f, nz, &ok);
		*out = (oi || oj || ok) ? 1 : 0;   // one flag for the triple, as the reference reports it
#endif
	}
	return (float)vol[i + (long long)j * nx + (long long)k * nx * ny];
}


// PERFORMANCE. This function is ~88% of the whole op's runtime (measured: 5318 of 6071
// sampling-profile hits), so it is written for the machine rather than for the eye. Three
// things the obvious transcription gets wrong, each measured separately:
//
//   1. `istep` was a LOOP-CARRIED dependency. Written as AFNI writes it --
//      `while (istep <= d1 && !out) { ...; if (!out) istep++; }` -- the next step's address
//      cannot be computed until this step's out-of-range flag is known, so every sample
//      pays the full int->float->fma->float->int->compare latency in series (~14 cycles of
//      pure chain) and nothing overlaps. The `if (!out) istep++` at the bottom plus `!out`
//      in the condition is EXACTLY `break on out` with a plain induction variable: when out
//      is set istep is not incremented and the condition fails immediately, so the loop ends
//      with istep unchanged either way. Rewritten as `for (istep...) { ...; if (out) break; }`
//      the addresses of all 21/16 steps are independent, the out-of-order engine overlaps
//      the gathers, and the only chain left is the min/max select. `out` is safe to test this
//      way because ss_sample_nn ASSIGNS it every call -- it never accumulates.
//   2. `o` may alias `p`, `n` and `st` as far as the compiler knows, so every store through
//      it (`o->means[1] +=`, `o->overshish[istep] =`, `o->n_over =`) forced a RELOAD of
//      p[0..2], n[0..2] and st->t on the next step -- seven loads per ray step of values
//      that never change. Everything the walk needs is pulled into locals first and only
//      the finished results are written through `o`.
//   3. The out-of-range test is one unsigned compare per axis, not a signed pair. Clamping
//      is done ONLY on the step that actually left the grid (that step is used and then the
//      walk stops), so the common in-range path carries no select chain at all.
//
// All three of those are bit-identical to the pre-optimisation release and are compiled into
// BOTH specialisations, verified on all nine benchmark images.
//
// The `memset` is the one thing that is NOT shared, and it is not there for correctness:
// nothing reads overshish past n_over, verified by poisoning the tail (0x7f fill -- all nine
// masks came back byte-identical). What it does is pin `o` in memory. Without it SROA promotes
// the struct and clang re-associates the surrounding float math enough to move the result
// (Dice 0.970-0.996), at ~4%. Faithful keeps it so it reproduces the release exactly; SS_FAST
// drops it. Note this one CANNOT be expressed as a runtime `if` even in principle: its job
// is to force `o`'s address to be taken so SROA cannot promote the struct, which is a
// compile-time property of the function, not a value you can branch on. Quality is the same either way -- against AFNI the two score 0.9407 and 0.9409
// mean Dice over the nine-image set. (Re-measure BOTH numbers on the SHIPPED pair whenever
// either kernel changes: 0.9413 was an intermediate build's figure and outlived it as a stale
// comment, then got copied into AGENTS.md as fact.)
static void SS_FN(ss_probe)(const unsigned char *vol, int nx, int ny, int nz,
		const ss_stats *st, const float *p, const float *n, int d1, int d4,
		ss_probe_t *o) {
	const int d2 = d1 / 2;   // AFNI: istep2max = ceil((d1/2)/travstp), travstp = 1 mm
	const float p0 = p[0], p1 = p[1], p2 = p[2];
	const float n0 = n[0], n1 = n[1], n2 = n[2];
	const float s_t = st->t, s_t2 = st->t2, s_tm = st->tm;
	const int nover_max = (d4 < 63) ? d4 : 63;
	float lmin = s_tm, lmax = s_t, lmind = 0.0f, lmaxd = 0.0f;
	float acc = 0.0f;
	int nm = 0, stopint = 0, out = 0, istep, nover = 0;

#ifndef SS_FAST
	memset(o, 0, sizeof(*o));   // pins `o` in memory; see the note above -- faithful mode only
#endif
	for (istep = 0; istep <= d1; istep++) {
		float val = SS_FN(ss_sample_nn)(vol, nx, ny, nz, p0 - n0 * istep, p1 - n1 * istep,
				p2 - n2 * istep, &out);
		if (val > s_t && !stopint) { acc += val; nm++; }
		else stopint = 1;
		if (lmin > val) { lmin = val; lmind = (float)istep; }
		if (istep <= d2 && lmax < val) { lmax = val; lmaxd = (float)istep; }
		if (out) break;
	}
	o->means[1] = acc / (float)nm;   // unguarded; convention 4 at ss_probe_t in skullstrip.c
	o->mm[0] = (s_t2 > lmin) ? s_t2 : lmin;
	o->mm[1] = (s_tm < lmax) ? s_tm : lmax;
	o->mmd[0] = lmind;
	o->mmd[1] = lmaxd;

	lmin = s_tm; lmax = s_t; lmind = lmaxd = 0.0f;
	acc = 0.0f; nm = 0; stopint = 0; out = 0;
#ifdef SS_FAST
	// Faithful's memset covers this; SS_FAST has to state it. Both call sites pass d4 = 15, so
	// the loop below always runs at least once and assigns means[0] -- but this header exists to
	// be included twice and invites reuse, and a caller passing d4 < 0 would read stack garbage
	// here (measured against a poisoned struct: -391.529 where faithful gives 0.0). means[0] is
	// the ONLY field at risk: `nover` is already 0-initialised so n_over is correct, and nothing
	// reads overshish past n_over.
	o->means[0] = 0.0f;
#endif
	for (istep = 0; istep <= nover_max; istep++) {
		float val = SS_FN(ss_sample_nn)(vol, nx, ny, nz, p0 + n0 * istep, p1 + n1 * istep,
				p2 + n2 * istep, &out);
		if (!istep) o->means[0] = val;
		o->overshish[istep] = val;
		nover = istep + 1;   // hoisted out of the store-through-o set; same final value
		if (val > s_t && !stopint) { acc += val; nm++; }
		else stopint = 1;
		if (lmin > val) { lmin = val; lmind = (float)istep; }
		if (lmax < val) { lmax = val; lmaxd = (float)istep; }
		if (out) break;
	}
	o->n_over = nover;
	o->means[2] = acc / (float)nm;
	o->over[0] = lmin; o->overd[0] = lmind;
	o->over[1] = lmax; o->overd[1] = lmaxd;
}

// it0..nit is a RANGE, not a count: AFNI re-enters this loop for each touchup stage
// with the SAME mesh and the SAME per-node ztv, advancing it0/nit each time. The
// temporal ramp lztfac is a function of absolute iteration, so continuing a stage is
// not the same as restarting one.
#ifdef SS_LINKAGE
int
#else
static int
#endif
SS_FN(ss_deform_range)(const unsigned char *vol, int nx, int ny, int nz, const ss_stats *st,
		ss_mesh *m, float *ztv, const float *stop, int it0, int nit, int niter_total,
		int nnsmooth, float *maxexp_out) {
	float *next = NULL, *smbuf = NULL;
	double l_mean;
	// BET's curvature sigmoid constants, in mm.
	const float rmin = 3.33f, rmax = 10.0f;
	const float E = 0.5f * (1.0f / rmin + 1.0f / rmax);
	const float F = 6.0f / (1.0f / rmin - 1.0f / rmax);
	const int d1 = (int)SS_D1_MM;
	const int d4 = 15;              // AFNI Opt->d4 = 15.0/rat, rat = 1 for human
	const float su1 = 1.0f;         // AFNI Opt->su1 = 1.0 -- NOT BET's 0.5
	// AFNI's smoothing schedule, and it is nothing like a per-N-iterations nudge:
	// NNsmooth = 72 geometric passes fired ONLY when it % smootheach == 0, it > 0 and
	// it < 0.75*N_it -- three events in the base 250-iteration stage (50, 100, 150),
	// and none at all once the stage loop pushes past 187. "Cannot filter too close to
	// convergence" is AFNI's own comment on the upper bound.
	// AFNI's -smootheach, -pushout and -avoid_eyes, all fixed at their AFNI values. These were
	// environment knobs while the constants were still being established; they are constants
	// now, because a stray SS_PUSHOUT=0 in a shell silently changes a shipped op's output with
	// no trace in that output -- exactly the reproducibility hazard the benchmark notes warn
	// about. SS_VAR and SS_NODE_DBG stay: the manifest's experiments depend on them.
	const int smootheach = 50;
	const int use_expansion = 1;
	const int no_eyes = 1;
	int var_lzt = ss_env_int("SS_VAR", 1, 0, 1);
	int node_dbg = ss_env_int("SS_NODE_DBG", -1, -1, 1 << 24);
	float threshclip;
	// FIXED for the whole run: AFNI's SO->Center is written once, when the icosphere is
	// created, and SUMA_RECOMPUTE_NORMALS never touches it. Tracking the running mesh
	// centroid instead makes both the inferior anti-leak floor and the eye zone drift
	// with the surface.
	const float *ctr = st->cog;
	double maxexp = 0.0;

	if (!vol || !st || !m || !ztv)
		return 1;
	next = (float *)malloc(sizeof(float) * 3 * (size_t)m->nv);
	smbuf = (float *)malloc(sizeof(float) * 6 * (size_t)m->nv);
	if (!next || !smbuf) { free(next); free(smbuf); return 1; }
	threshclip = (st->t98 - st->t2) / 5.0f;
	ss_mesh_normals(m);

	for (int it = it0; it < nit; it++) {
		float lztfac;
		// Expand fast early to escape large CSF pools near the outer surface, then
		// tighten towards the end. Ramps 0.2 -> 1.2 over the BASE iteration budget and
		// pins at 1.2 once the stage loop runs past it.
		if (var_lzt) {
			if (it <= niter_total) {
				lztfac = 1.2f * (float)it / (float)niter_total;
				if (lztfac < 0.2f) lztfac = 0.2f;
			} else
				lztfac = 1.2f;
		} else
			lztfac = 1.0f;

		l_mean = ss_mean_seg_len(m);
		maxexp = 0.0;

#ifdef SS_FAST
		// Every node is independent: it reads the PREVIOUS iteration's m->v/m->nrm plus
		// read-only ztv/stop/vol/st and writes only its own next[3i..3i+2]. The single
		// cross-node value, maxexp, is a MAX reduction -- exact and order-independent on
		// doubles -- so this is bit-identical at any thread count, which the benchmark set is
		// checked against. It is NOT a sum: a sum reduction here would make output
		// thread-dependent. Serialised when the SS_NODE_DBG dev hook is armed so its per-node
		// trace stays in node order and stays diffable against AFNI's dbg_<node>.1D.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max : maxexp) if (node_dbg < 0)
#endif
#endif
		for (int i = 0; i < m->nv; i++) {
			const float *p = m->v + 3 * i, *n = m->nrm + 3 * i;
			float sx = 0, sy = 0, sz = 0;
			float sv[3], sn[3], stg[3], dp, nsn;
			float su2, su3, su4 = 0.0f, lZt, tb, f3, r;
			int deg = m->nbr_off[i + 1] - m->nbr_off[i];
			float Imin, Imax;
			ss_probe_t pr;
			float max_normal_exp = -1.0f;   // <0 means "no cap"
			double un;

			if (deg <= 0) {
				memcpy(next + 3 * i, p, sizeof(float) * 3);
				continue;
			}
			for (int q = m->nbr_off[i]; q < m->nbr_off[i + 1]; q++) {
				const float *pj = m->v + 3 * m->nbr[q];
				sx += pj[0]; sy += pj[1]; sz += pj[2];
			}
			sv[0] = sx / deg - p[0]; sv[1] = sy / deg - p[1]; sv[2] = sz / deg - p[2];
			dp = sv[0] * n[0] + sv[1] * n[1] + sv[2] * n[2];
			for (int a = 0; a < 3; a++) {
				sn[a] = dp * n[a];
				stg[a] = sv[a] - sn[a];
			}
			nsn = sqrtf(sn[0] * sn[0] + sn[1] * sn[1] + sn[2] * sn[2]);
			r = (nsn > 1e-8f) ? (float)(l_mean * l_mean) / (2.0f * nsn) : 1e8f;
			su2 = 0.5f * (1.0f + tanhf(F * (1.0f / r - E)));

			SS_FN(ss_probe)(vol, nx, ny, nz, st, p, n, d1, d4, &pr);
			Imin = pr.mm[0];
			Imax = pr.mm[1];

			// Local shrink factor: per-node base times the temporal ramp, with a
			// SPATIAL floor for nodes well below the surface centre. That inferior
			// clip is what stops the surface running down into the neck, and it is a
			// position test -- not a schedule.
			if (p[2] - ctr[2] < -25.0f)
				lZt = (SS_SHRINK_BOT_NOEDGE > ztv[i] * lztfac) ? SS_SHRINK_BOT_NOEDGE
				                                               : ztv[i] * lztfac;
			else
				lZt = ztv[i] * lztfac;
			if (lZt > 1.0f)
				lZt = 0.999f;

			tb = (Imax - st->t2) * lZt + st->t2;
			f3 = 2.0f * (Imin - tb) / (Imax - st->t2);
			su3 = SS_EXP_FRAC * (float)l_mean * f3;

			if (use_expansion) {
				float tbe = (pr.over[1] - st->t2) * lZt + st->t2;
				float f4 = 2.0f * (pr.over[0] - tbe) / (pr.over[1] - st->t2);
				su4 = SS_EXP_FRAC * (float)l_mean * f4;
				if (su4 < 0.0f)
					su4 = 0.0f; // outward-only: it expands, never pulls back
			}

			// Eye avoidance. In the anterior-inferior octant relative to the surface
			// centre, a node whose mean ABOVE clearly exceeds its mean BELOW is
			// looking into the eye, so cap how far it may expand -- at the first
			// significant jump in the outward profile, which is the eye's leading edge.
			if (no_eyes && it > (int)(0.1f * niter_total) &&
					SS_IS_EYE_ZONE(p, ctr) && !(stop && stop[i] != 0.0f)) {
				if (pr.means[2] > 1.2f * pr.means[1]) {
					max_normal_exp =
						(float)ss_max_shish_jump(pr.overshish, pr.n_over, threshclip);
				} else if (pr.means[0] >= st->t98) {
					max_normal_exp = 0.0f;
				}
			}
			// A frozen node (Stop < 0) stops expanding entirely for the rest of the
			// run. This is AFNI's brake on leakage into fat/skull, and it is the half
			// of the touchup feedback that a relax-only implementation misses.
			if (stop && stop[i] < 0.0f) {
				su3 = 0.0f;
				su4 = 0.0f;
			}
			// Cap total normal growth. This is the brake that stops a blurry volume
			// from ballooning: without it the outward-only su4 has nothing opposing it.
			if (max_normal_exp >= 0.0f && (su3 + su4) > max_normal_exp) {
				su3 = max_normal_exp;
				su4 = 0.0f;
			}

			if (node_dbg == i) {
				// Column-for-column diffable against AFNI's dbg_<node>.1D (-node_dbg N
				// -debug 3). Same order as its header block.
				double dbg[38] = {
					l_mean, Imin, Imax, tb,
					pr.mmd[0], pr.mmd[1], pr.over[0], pr.over[1],
					pr.overd[0], pr.overd[1],
					pr.means[0], pr.means[1], pr.means[2],
					p[0], p[1], p[2], dp,
					sx / deg, sy / deg, sz / deg,
					n[0], n[1], n[2],
					sv[0], sv[1], sv[2],
					stg[0], stg[1], stg[2],
					sn[0], sn[1], sn[2],
					su1, su2, su3,
					su1 * stg[0] + su2 * sn[0] + (su3 + su4) * n[0],
					su1 * stg[1] + su2 * sn[1] + (su3 + su4) * n[1],
					su1 * stg[2] + su2 * sn[2] + (su3 + su4) * n[2]};
				fprintf(stderr, "%d %d", i, it);
				for (int q = 0; q < 38; q++)
					fprintf(stderr, " %.6f", dbg[q]);
				fprintf(stderr, " %.6f %.6f %.6f %.6f\n",
						(double)st->t2, (double)st->t, (double)st->tm, (double)st->t98);
			}

			for (int a = 0; a < 3; a++)
				next[3 * i + a] = p[a] + su1 * stg[a] + su2 * sn[a] + (su3 + su4) * n[a];
			// Never let a node leave the working grid (AFNI's SUMA_LimitCoordToVolume;
			// it fires on unusual anatomy, and without it the ray walk starts clamping).
			// Written so NaN is caught, not passed through: `x < 0` is FALSE for NaN, so a
			// plain pair of comparisons lets a NaN coordinate survive into the next
			// iteration's ss_sample_nn, where (int)(NaN + 0.49f) is undefined (arm64 gives
			// 0, wasm32 traps). `!(x >= 0)` is true for NaN.
			for (int a = 0; a < 3; a++) {
				float hi = (a == 0) ? nx - 1.0f : ((a == 1) ? ny - 1.0f : nz - 1.0f);
				if (!(next[3 * i + a] >= 0.0f)) next[3 * i + a] = 0.0f;
				else if (next[3 * i + a] > hi) next[3 * i + a] = hi;
			}
			un = (double)(su1 * stg[0] + su2 * sn[0] + (su3 + su4) * n[0]);
			{
				double uy = (double)(su1 * stg[1] + su2 * sn[1] + (su3 + su4) * n[1]);
				double uz = (double)(su1 * stg[2] + su2 * sn[2] + (su3 + su4) * n[2]);
				double mag = sqrt(un * un + uy * uy + uz * uz);
				if (mag > maxexp) maxexp = mag;
			}
		}
		memcpy(m->v, next, sizeof(float) * 3 * (size_t)m->nv);

		if (nnsmooth > 0 && smootheach > 0 && (it % smootheach) == 0 && it > 0 &&
				it < 0.75f * niter_total)
			ss_smooth_nn_wrap(m, nnsmooth, ctr, smbuf);

		ss_mesh_normals(m);
	}
	if (maxexp_out) *maxexp_out = (float)maxexp;
	free(next);
	free(smbuf);
	return 0;
}

