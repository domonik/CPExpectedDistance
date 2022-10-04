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


static void reprint(PyObject *obj) {
    PyObject* repr = PyObject_Repr(obj);
    PyObject* str = PyUnicode_AsEncodedString(repr, "utf-8", "~E~");
    const char *bytes = PyBytes_AS_STRING(str);

    printf("REPR: %s\n", bytes);

    Py_XDECREF(repr);
    Py_XDECREF(str);
}


static int set_config_double(PyObject *md_config, double *value, char *key){
    PyObject *cur_value = PyDict_GetItemString(md_config, key);
    if (cur_value){
        if (PyFloat_Check(cur_value)){
        }
        else if (PyLong_Check(cur_value)) {
            cur_value = PyNumber_Float(cur_value);
            if (!cur_value){
                PyErr_SetString(PyExc_TypeError, "Expected Float got something else");
                return 0;
            }
        }
        else{
            PyErr_SetString(PyExc_TypeError, "Expected Float got something else");
            return 0;
        }
        double v = PyFloat_AS_DOUBLE(cur_value);
        *value = v;
    }
    Py_XDECREF(cur_value);
    return 1;
}

static int set_config_int(PyObject *md_config, int *value, char *key){
    PyObject *cur_value = PyDict_GetItemString(md_config, key);
    if (cur_value){
        if (!PyLong_Check(cur_value)){
            PyErr_SetString(PyExc_TypeError, "Expected Integer got something else");
            return 0;
        }
        int v = (int)PyLong_AsLong(cur_value);
        *value = v;
    }
    Py_XDECREF(cur_value);
    return 1;
}


static int set_model_details_from_config(PyObject *md_config, vrna_md_t *mdt){
    if (!set_config_double(md_config, &mdt->temperature, "temperature")) {return 0;};
    if (!set_config_int(md_config, &mdt->dangles, "dangles")) {return 0;};
    if (!set_config_int(md_config, &mdt->pf_smooth, "pf_smooth")) {return 0;};
    if (!set_config_int(md_config, &mdt->special_hp, "special_hp")) {return 0;};
    if (!set_config_int(md_config, &mdt->noLP, "noLP")) {return 0;};
    if (!set_config_int(md_config, &mdt->noGU, "noGU")) {return 0;};
    if (!set_config_int(md_config, &mdt->noGUclosure, "noGUclosure")) {return 0;};
    if (!set_config_int(md_config, &mdt->logML, "logML")) {return 0;};
    if (!set_config_int(md_config, &mdt->circ, "circ")) {return 0;};
    if (!set_config_int(md_config, &mdt->gquad, "gquad")) {return 0;};
    if (!set_config_int(md_config, &mdt->uniq_ML, "uniq_ML")) {return 0;};
    if (!set_config_int(md_config, &mdt->energy_set, "energy_set")) {return 0;};
    if (!set_config_int(md_config, &mdt->backtrack, "backtrack")) {return 0;};
    if (!set_config_int(md_config, &mdt->compute_bpp, "compute_bpp")) {return 0;};
    if (!set_config_int(md_config, &mdt->max_bp_span, "max_bp_span")) {return 0;};
    if (!set_config_int(md_config, &mdt->min_loop_size, "min_loop_size")) {return 0;};
    if (!set_config_int(md_config, &mdt->window_size, "window_size")) {return 0;};
    if (!set_config_int(md_config, &mdt->oldAliEn, "oldAliEn")) {return 0;};
    if (!set_config_int(md_config, &mdt->ribo, "ribo")) {return 0;};
    if (!set_config_double(md_config, &mdt->betaScale, "betaScale")) {return 0;};
    if (!set_config_double(md_config, &mdt->cv_fact, "cv_fact")) {return 0;};
    if (!set_config_double(md_config, &mdt->nc_fact, "nc_fact")) {return 0;};
    if (!set_config_double(md_config, &mdt->sfact, "sfact")) {return 0;};
    return 1;
}


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



static PyObject *clote_ponty_method_old(PyObject *self, PyObject *args) {
    PyObject *md_config = NULL;
    double z;
    double *res;
    double *prev_res;
    char *seq = NULL;
    vrna_md_t md;
    vrna_md_set_default(&md);



    /* Parse arguments */

    if(!PyArg_ParseTuple(args, "sO", &seq, &md_config)) {

        return NULL;

    }
    if (!PyDict_Check(md_config)){
        PyErr_SetString(PyExc_TypeError, "Second Argument needs to be a Dictionary");
        return NULL;
    }
//    printf("temp: %f\n", md.temperature);
    if (!set_model_details_from_config(md_config, &md)) {return NULL;};
//    printf("temp: %d\n", md.dangles);


    vrna_fold_compound_t  *fc = vrna_fold_compound(seq, &md, VRNA_OPTION_DEFAULT);
    double mfe = (double)vrna_mfe(fc, NULL);
    vrna_exp_params_rescale(fc, &mfe);

    vrna_pf(fc, NULL);
    long int dims[2] = {fc->exp_matrices->length, fc->exp_matrices->length };
    PyArrayObject *output;
    output = (PyArrayObject *) PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);

    unsigned int j;
    for (unsigned int l = 1; l < fc->exp_matrices->length + 1; l = l+1){
        for (unsigned int i = 1;  i < fc->exp_matrices->length + 1 - l; i = i+1){
            j = i + l;
            res = PyArray_GETPTR2(output, (i-1), (j-1));
            prev_res =  PyArray_GETPTR2(output, (i-1), (j-2));
            z = *prev_res * fc->exp_matrices->scale[1] + get_Z(fc, i, j);
            /*if (i == 1 && j == 10){
                printf("i: %d j: %d z: %f\n", i, j,  z);
                printf("prev: %f sc: %f z: %f\n", *prev_res, (double)fc->exp_matrices->scale[1],  get_Z(fc, i, j));
            }*/
            for (unsigned int  k = i+1;  k <= j; k = k+1){
                prev_res = PyArray_GETPTR2(output, i-1, k-2);
                z +=  (*prev_res + get_Z(fc, i, k-1)) * get_ZB(fc, k, j) * new_exp_E_ext_stem(fc, k, j);
                /*if (i == 1 && j == 10){
                    printf("k: %d", k);
                    printf("prev_i: %f sc: %f z: %f exp: %f\n", *prev_res, get_Z(fc, i, k-1),  get_ZB(fc, k, j), new_exp_E_ext_stem(fc, k, j));
                }*/

            }
//                        printf("i: %d j: %d z: %f\n", i, j, z);

            *res = z ;

        }
    }
    for (unsigned int i = 1; i < fc->exp_matrices->length + 1; i = i+1){
        for (unsigned int j = i+1;  j < fc->exp_matrices->length + 1; j = j+1){
            res = PyArray_GETPTR2(output, i-1, j-1);
            *res /= get_Z(fc, i, j);
        }
    }
    vrna_fold_compound_free(fc);

    return PyArray_Return(output);
}

static PyMethodDef ExpDMethods[] = {
    {"cp_expected_distance_old", clote_ponty_method_old, METH_VARARGS, "Python interface for clote-ponty expected distance calculation"},
    {"cp_expected_distance", clote_ponty_method, METH_VARARGS, "Python interface for getting the fold compound back"},

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
