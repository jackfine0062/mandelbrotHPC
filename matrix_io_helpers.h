#include <string.h>
#include <errno.h>
#include <ctype.h>

////////// CSV File Reading //////////

static inline float __read_csv_val(const char* tok) {
    float val;
    int off;
    if (sscanf(tok, " %f %n", &val, &off) != 1 || off != strlen(tok)) {
        fprintf(stderr, "Not a number in CSV file, using 0.0: %s\n", tok);
        return 0.0;
    }
    return val;
}

static inline size_t __read_csv_first_line(char* line, float** out) {
    char* tok = strtok(line, ",\n");
    size_t i = 0, size = 16;
    float* vals = (float*)malloc(16*sizeof(float));
    while (tok) {
        vals[i] = __read_csv_val(tok);
        if (++i == size) {
            size += size > 512 ? 512 : size;
            vals = (float*)realloc(vals, size*sizeof(float));
        }
        tok = strtok(NULL, ",");
    }
    *out = vals;
    return i;
}

static inline void __read_csv_line(char* line, float* out, size_t count) {
    char* tok = strtok(line, ",\n");
    size_t i = 0;
    while (tok) {
        out[i] = __read_csv_val(tok);
        if (++i == count) { return; } // stop reading
        tok = strtok(NULL, ",");
    }
    memset(out+i, 0, (count-i)*sizeof(float)); // zero-fill remainder
}


////////// NPY File Reading //////////

static inline const char* __py_dict_value(const char* dict, const char* key) {
    size_t key_len = strlen(key);
    key = strstr(dict, key);
    if (key == NULL || key[-1] != key[key_len] ||
        (key[key_len] != '\'' && key[key_len] != '"')) {
        return NULL;
    }
    const char* s = key + key_len + 1;
    while (isspace(*s)) { s++; }
    if (*s != ':') { return NULL; }
    s++;
    while (isspace(*s)) { s++; }
    return s;
}

static inline char* __py_dict_value_str(const char* dict, const char* key) {
    const char* s = __py_dict_value(dict, key);
    if (!s) { return NULL; }
    char c = *s;
    if (c != '\'' && c != '"') { return NULL; }
    const char* end = ++s;
    while (*end != c) { end++; }
    size_t len = end - s;
    char* data = malloc(len+1);
    strncpy(data, s, len);
    data[len] = 0;
    return data;
}

static inline bool __py_dict_value_bool(const char* dict,
                                        const char* key, bool* val) {
    const char* s = __py_dict_value(dict, key);
    if (s) {
        if (strncmp(s, "False", 5) == 0) { *val = false; return true; }
        else if (strncmp(s, "True", 4) == 0) { *val = true; return true; }
    }
    return false;
}

static inline bool __py_dict_value_tuple(const char* dict,
                                         const char* key, size_t* val) {
    const char* s = __py_dict_value(dict, key);
    val[0] = val[1] = 1;
    char c = 0;
    return s && *s++ == '(' && ((sscanf(s, " %c", &c) == 1 && c == ')') ||
        (sscanf(s, " %zu %c", &val[0], &c) == 2 && c == ')') ||
        (sscanf(s, " %zu , %c", &val[0], &c) == 2 && c == ')') ||
        (sscanf(s, " %zu , %zu %c", &val[0], &val[1], &c) == 3 && c == ')') ||
        (sscanf(s, " %zu , %zu , %c", &val[0], &val[1], &c) == 3 && c == ')'));
}

static inline bool __npy_read_header(FILE* file, size_t* sh, size_t* offset) {
    unsigned char header[10];
    if (fread(header, 1, 10, file) != 10) { return false; }
    if (memcmp(header, "\x93NUMPY", 6) != 0) { errno = EINVAL; return false; }
    // header[6] is major file version
    // header[7] is minor file version
    int len = *(unsigned short*)(header+8); // assumes running on little-endian
    *offset = sizeof(header) + len;
    char* dict = (char*)malloc(len+1);
    if (fread(dict, 1, len, file) != len || dict[0] != '{') {
        free(dict);
        errno = EINVAL;
        return false;
    }
    dict[len] = 0;

    // only allowed descr: 'float32', 'f4', or '<f4'
    char* descr = __py_dict_value_str(dict, "descr");
    if (!descr) { free(dict); errno = EINVAL; return false; }
    if (strcmp(descr, "float32") != 0 && strcmp(descr, "f4") != 0 &&
        strcmp(descr, "<f4") != 0) {
        errno = EINVAL;
        free(descr);
        free(dict);
        return false;
    }
    free(descr);

    // only allowed fortran_order is false
    bool fortran_order = false;
    if (!__py_dict_value_bool(dict, "fortran_order", &fortran_order) &&
        fortran_order) {
        errno = EINVAL;
        free(dict);
        return false;
    }

    // only allowed to be 0d, 1d, or 2d, but this is checked elsewhere
    if (!__py_dict_value_tuple(dict, "shape", sh) || sh[0] < 1 || sh[1] < 1) {
        errno = EINVAL;
        free(dict);
        return false;
    }
    free(dict);
    return true;
}
