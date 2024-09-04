#include <string.h>

// Define a default GIT_TAG if not provided
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "Unknown"
#endif

extern "C" {
    void get_version_tag_c(char* version, int* length)
    {
        const char* tag = PACKAGE_VERSION;
        *length = strlen(tag);
        strcpy(version, tag);
    }
}