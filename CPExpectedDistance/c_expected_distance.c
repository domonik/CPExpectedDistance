#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <stdlib.h>
#include <stdio.h>

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/model.h>
#include <numpy/arrayobject.h>
#include <ViennaRNA/loops/external.h>
#include <ViennaRNA/alphabet.h>

typedef struct {
    PyObject_HEAD
    void *ptr;
    void *ty;
    int own;
    PyObject *next;
    PyObject *dict;
} SwigPyObject;



static double get_Z(vrna_fold_compound_t *fc, int i, int j){
    int idx = (((fc->exp_matrices->length + 1 - i) * (fc->exp_matrices->length - i)) / 2) + fc->exp_matrices->length + 1 - j;
    return fc->exp_matrices->q[idx];
}

static double get_ZB(vrna_fold_compound_t *fc, int i, int j){
    int idx = (((fc->exp_matrices->length + 1 - i) * (fc->exp_matrices->length - i)) / 2) + fc->exp_matrices->length + 1 - j;
    return fc->exp_matrices->qb[idx];
}

static double new_exp_E_ext_stem(vrna_fold_compound_t *fc, unsigned int i, unsigned int j){
    unsigned int type;
    int enc5, enc3;
    enc5 = enc3 = -1;

    type = vrna_get_ptype_md(fc->sequence_encoding2[i],
                             fc->sequence_encoding2[j],
                             &(fc->params->model_details));

    if (i > 1)
      enc5 = fc->sequence_encoding[i - 1];
    if (j < fc->length)
      enc3 = fc->sequence_encoding[j + 1];

    return (double)vrna_exp_E_ext_stem(type,enc5,enc3,fc->exp_params);
}

static vrna_fold_compound_t *swig_fc_to_fc(PyObject *swig_fold_compound) {
    SwigPyObject *swf = (SwigPyObject*) swig_fold_compound;
    vrna_fold_compound_t *fc = (vrna_fold_compound_t*) swf->ptr;
    return fc;
}

static void fill_expected_distance(PyArrayObject *parray, vrna_fold_compound_t *fc) {
    double *res;
    double z;
    unsigned int j;
    double *prev_res;

    for (unsigned int l = 1; l < fc->exp_matrices->length + 1; l = l+1){
        for (unsigned int i = 1;  i < fc->exp_matrices->length + 1 - l; i = i+1){
            j = i + l;
            res = PyArray_GETPTR2(parray, (i-1), (j-1));
            prev_res =  PyArray_GETPTR2(parray, (i-1), (j-2));
            z = *prev_res * fc->exp_matrices->scale[1] + get_Z(fc, i, j);
            for (unsigned int  k = i+1;  k <= j; k = k+1){
                prev_res = PyArray_GETPTR2(parray, i-1, k-2);
                z +=  (*prev_res + get_Z(fc, i, k-1)) * get_ZB(fc, k, j) * new_exp_E_ext_stem(fc, k, j);
            }

            *res = z ;
        }
    }
    for (unsigned int i = 1; i < fc->exp_matrices->length + 1; i = i+1){
        for (unsigned int j = i+1;  j < fc->exp_matrices->length + 1; j = j+1){
            res = PyArray_GETPTR2(parray, i-1, j-1);
            *res /= get_Z(fc, i, j);
        }
    }
}



static PyObject *clote_ponty_method(PyObject *self, PyObject *args) {
    PyObject *swig_fc;

    if(!PyArg_ParseTuple(args, "O", &swig_fc)) {

    return NULL;

    }
    vrna_fold_compound_t *fc = swig_fc_to_fc(swig_fc);
    double mfe = (double)vrna_mfe(fc, NULL);
    vrna_exp_params_rescale(fc, &mfe);
    vrna_pf(fc, NULL);
    long int dims[2] = {fc->exp_matrices->length, fc->exp_matrices->length };
    PyArrayObject *output = (PyArrayObject *) PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);
    fill_expected_distance(output, fc);
    return PyArray_Return(output);
}


static PyMethodDef ExpDMethods[] = {
    {"cp_expected_distance", clote_ponty_method, METH_VARARGS, "Python interface for computing clote-ponty expected distance"},

    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef c_expected_distance = {
    PyModuleDef_HEAD_INIT,
    "c_expected_distance",
    "clote-ponty expected distance calculation",
    -1,
    ExpDMethods
};

PyMODINIT_FUNC PyInit_c_expected_distance(void) {
    import_array();
    return PyModule_Create(&c_expected_distance);
}
