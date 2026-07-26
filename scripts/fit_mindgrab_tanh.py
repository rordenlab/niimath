#!/usr/bin/env python3
"""Derive the float32 tanh rational approximation used by src/mindgrab.c.

MindGrab evaluates tanh-GELU 6.3e9 times per volume, so libm tanhf dominates the
runtime. This fits an odd rational  x * P(x^2) / Q(x^2)  on [-9, 9] (tanhf rounds to
+-1 beyond that) by the usual linearised rational least squares -- solve for
x*P(u) - t*Q(u) = 0 weighted by 1/Q from the previous iterate -- then reports the
worst float32 ULP error against libm.

Run it to reproduce or re-tune MG_TANH_* in mindgrab.c; it is not part of any build.
"""
import numpy as np

LIM = 9.0
NP, NQ = 7, 4  # coefficients of P and Q in powers of u = x^2 (Q0 pinned to 1)


def fit(iters=400, boost=0.35):
    # dense sampling, denser near 0 where the relative error matters most
    x = np.concatenate([np.linspace(1e-9, 1.0, 60000), np.linspace(1.0, LIM, 140000)])
    t = np.tanh(x)
    # fit in the normalised variable m = (x/LIM)^2 in [0,1]; u^6 = 2.8e11 otherwise
    # and the design matrix is too ill-conditioned to solve.
    u = (x / LIM) ** 2
    A = np.empty((x.size, NP + NQ - 1))
    for i in range(NP):
        A[:, i] = x * u ** i
    for j in range(1, NQ):
        A[:, NP + j - 1] = -t * u ** j
    qprev = np.ones_like(x)
    w = np.ones_like(x)
    best = None
    for it in range(iters):
        s = w / qprev
        sol, *_ = np.linalg.lstsq(A * s[:, None], (t * s), rcond=None)
        p = sol[:NP]
        q = np.concatenate([[1.0], sol[NP:]])
        qv = np.polyval(q[::-1], u)
        if np.any(qv <= 1e-6):
            break  # denominator lost positivity; keep the last good iterate
        qprev = qv
        err = np.abs(x * np.polyval(p[::-1], u) / qv - t)
        if best is None or err.max() < best[0]:
            best = (err.max(), p, q)
        # gentle Remez-style equalisation of the absolute error
        w = w * (1.0 + boost * err / err.max())
        w /= w.max()
    if best is None:
        raise SystemExit("fit diverged")
    print("# fitted max abs err (float64) %.3e" % best[0])
    scale = np.array([LIM ** (-2 * i) for i in range(max(NP, NQ))])
    return best[1] * scale[:NP], best[2] * scale[:NQ]


def report(p, q):
    p32 = np.float32(p)
    q32 = np.float32(q)

    def approx(x):
        x = np.clip(np.float32(x), np.float32(-LIM), np.float32(LIM))
        u = np.float32(x * x)
        pr = np.float32(np.float32(p32[-1]) + np.zeros_like(u))
        for c in p32[-2::-1]:
            pr = np.float32(pr * u + c)
        qr = np.float32(np.float32(q32[-1]) + np.zeros_like(u))
        for c in q32[-2::-1]:
            qr = np.float32(qr * u + c)
        return np.float32(np.float32(x * pr) / qr)

    xs = np.float32(np.unique(np.concatenate([
        np.linspace(-12, 12, 600001), np.linspace(-1.5, 1.5, 300001)])))
    ref = np.tanh(xs.astype(np.float64)).astype(np.float32)
    got = approx(xs)
    ai = ref.view(np.int32).astype(np.int64)
    bi = got.view(np.int32).astype(np.int64)
    ai = np.where(ai < 0, np.int64(-2147483648) - ai, ai)
    bi = np.where(bi < 0, np.int64(-2147483648) - bi, bi)
    ulp = np.abs(ai - bi)
    print("max ulp %d   max abs err %.3e   worst at x = %g"
          % (ulp.max(), np.abs(got - ref).max(), xs[np.argmax(ulp)]))
    for n, c in (("P", p32), ("Q", q32)):
        for i, v in enumerate(c):
            print("#define MG_TANH_%s%d %.9gf" % (n, 2 * i + (1 if n == "P" else 0), v))


if __name__ == "__main__":
    report(*fit())
