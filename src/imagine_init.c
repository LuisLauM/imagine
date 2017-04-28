#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP imagine_engine1(SEXP, SEXP);
extern SEXP imagine_engine2(SEXP, SEXP);
extern SEXP imagine_engine3(SEXP, SEXP, SEXP, SEXP);
extern SEXP imagine_engine4(SEXP, SEXP);
extern SEXP imagine_engine5(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"imagine_engine1", (DL_FUNC) &imagine_engine1, 2},
    {"imagine_engine2", (DL_FUNC) &imagine_engine2, 2},
    {"imagine_engine3", (DL_FUNC) &imagine_engine3, 4},
    {"imagine_engine4", (DL_FUNC) &imagine_engine4, 2},
    {"imagine_engine5", (DL_FUNC) &imagine_engine5, 4},
    {NULL, NULL, 0}
};

void R_init_imagine(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
