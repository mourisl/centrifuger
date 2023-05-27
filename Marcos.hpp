#ifndef _MOURISL_MARCOS
#define _MOURISL_MARCOS

#define SAVE_VAR(fp, x) (fwrite(&(x), sizeof(x), 1, fp))
#define LOAD_VAR(fp, x) (fread(&(x), sizeof(x), 1, fp))

#endif
