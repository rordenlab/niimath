#ifndef AL_THREAD_LOCAL_H
#define AL_THREAD_LOCAL_H

/*
 * Thread-local storage used by the OpenMP callback/workspace state.
 *
 * Keep the compiler spellings here so shared sources do not each grow their own
 * portability assumptions. MSVC's C compiler does not accept GNU __thread;
 * C11 compilers use the standard spelling, with __thread retained as a fallback
 * for older GCC/Clang modes.
 *
 * Serial builds deliberately use ordinary static storage: no concurrent worker
 * can enter these callback/workspace paths when OpenMP is disabled.
 */
#if defined(_OPENMP)
#  if defined(_MSC_VER)
#    define AL_THREAD_LOCAL __declspec(thread)
#  elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#    define AL_THREAD_LOCAL _Thread_local
#  elif defined(__GNUC__) || defined(__clang__)
#    define AL_THREAD_LOCAL __thread
#  else
#    error "OpenMP build requires a supported thread-local storage qualifier"
#  endif
#else
#  define AL_THREAD_LOCAL
#endif

#endif /* AL_THREAD_LOCAL_H */
