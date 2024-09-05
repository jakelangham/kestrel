// This file is part of the Kestrel software for simulations
// of sediment-laden Earth surface flows.
//
// Version v1.1.1
//
// Copyright 2023 Mark J. Woodhouse, Jake Langham, (University of Bristol).
//
// This program is free software: you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free 
// Software Foundation, either version 3 of the License, or (at your option) 
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT 
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
// more details.
//
// You should have received a copy of the GNU General Public License along with 
// this program. If not, see <https://www.gnu.org/licenses/>. 

#include <string.h>

// Define a default GIT_TAG if not provided
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "v1.1.1"
#endif

extern "C" {
    void get_version_tag_c(char* version, int* length)
    {
        const char* tag = PACKAGE_VERSION;
        *length = strlen(tag);
        strcpy(version, tag);
    }
}
