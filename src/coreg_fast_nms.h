#ifndef COREG_FAST_NMS_H
#define COREG_FAST_NMS_H

#include <float.h>
#include <math.h>
#include <stddef.h>

/* Internal coarse-grid candidate representation. `order` is the deterministic
 * enumeration order and makes equal-cost sorting independent of qsort details. */
typedef struct {
    double p[12];
    double c;
    size_t order;
} cf_nms_candidate;

/* Historical rigid coarse paths keep a plain cost-sorted top-N without angular
 * suppression. Strict cost comparison preserves enumeration order on ties. */
static inline void cf_insert_sorted_candidate(cf_nms_candidate *top, int ntop,
                                              const double p[12], double cost,
                                              size_t order) {
    if (!top || !p || ntop <= 0) return;
    for (int t = 0; t < ntop; t++) {
        if (cost < top[t].c) {
            for (int u = ntop - 1; u > t; u--) top[u] = top[u - 1];
            for (int i = 0; i < 12; i++) top[t].p[i] = p[i];
            top[t].c = cost;
            top[t].order = order;
            break;
        }
    }
}

/*
 * Select the lowest-cost angularly distinct candidates.
 *
 * Selecting globally by cost before accepting each basin is important: replacing
 * an already-retained basin in place can break both cost ordering and the
 * minimum-separation invariant. Repeated scans produce the same greedy total order
 * without comparator-based qsort (notably expensive in WASM). At most six basins
 * are retained from the <=750-point grid, so this remains small and bounded.
 *
 * Equal-cost ties retain grid-enumeration order. Returns the number written to
 * `selected`; `pool` is unchanged.
 */
static inline int cf_select_angular_nms(const cf_nms_candidate *pool, size_t npool,
                                        cf_nms_candidate *selected, int max_selected,
                                        double min_deg, double param_to_deg,
                                        double penalty) {
    if (!pool || !selected || npool == 0 || max_selected <= 0 ||
        !(min_deg > 0.0 && min_deg <= DBL_MAX) ||
        !(param_to_deg > 0.0 && param_to_deg <= DBL_MAX))
        return 0;
    int nselected = 0;
    while (nselected < max_selected) {
        size_t best = npool;
        for (size_t i = 0; i < npool; i++) {
            double cv = pool[i].c;
            if (!(cv >= -DBL_MAX && cv <= DBL_MAX) || !(cv < penalty)) continue;
            int distinct = 1;
            for (int j = 0; j < nselected; j++) {
                double d0 = (pool[i].p[3] - selected[j].p[3]) * param_to_deg;
                double d1 = (pool[i].p[4] - selected[j].p[4]) * param_to_deg;
                double d2 = (pool[i].p[5] - selected[j].p[5]) * param_to_deg;
                if (sqrt(d0*d0 + d1*d1 + d2*d2) < min_deg) {
                    distinct = 0;
                    break;
                }
            }
            if (!distinct) continue;
            if (best == npool || cv < pool[best].c ||
                (cv == pool[best].c && pool[i].order < pool[best].order))
                best = i;
        }
        if (best == npool) break;
        selected[nselected++] = pool[best];
    }
    return nselected;
}

#endif /* COREG_FAST_NMS_H */
