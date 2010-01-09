/**
 *    AT_FileOperations.h
 *    ===================
 *
 *    Created on: 28.07.2009
 *    Author: greilich
 *
 *    Copyright 2006, 2009 Steffen Greilich / the libamtrack team
 *
 *    This file is part of the AmTrack program (libamtrack.sourceforge.net).
 *
 *    AmTrack is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    AmTrack is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with AmTrack (file: copying.txt).
 *    If not, see <http://www.gnu.org/licenses/>
 */

#include "AT_FileOperations.h"

void AT_browseInput(  char*  fileName,
            long*  nLines,
            long*  n_gamma_parameter)
{
  FILE*  inputFile;
  inputFile  =  fopen(  fileName, "r");
  char* retval;
  if (inputFile == 0) return;

    char line[256];

    *nLines    =  -1;
    while(fgets(line, 256, inputFile))
    {
    if(strstr(line, ":start spectrum:") != NULL)
    {
      while(strstr(line, ":stop spectrum:") == NULL)
      {
        retval = fgets(line, 256, inputFile);
        (*nLines)++;
      }
    }

    if(strstr(line, "model") != NULL)
    {
      char  seps[]  = "=";                  // delimiters can be ",", ";", or tab
      char*  token;

      token  =  strtok( line, seps );

      token    =  strtok( NULL, seps );
      long model  =  atol(token);

      if(model == 1){
        *n_gamma_parameter = 2;
      }
    }
  }
  return;
}

void AT_browseInputS(  char**  fileName,
            // return:
            long*  nLines,
            long*  n_gamma_parameter){
  AT_browseInput(  *fileName,
            // return:
            nLines,
            n_gamma_parameter);};

void AT_readInput(  char*  fileName,
          long*  nLines,
          long*  n_gamma_parameter,
          float*  E_MeV_u,
          long*  particle_no,
          float*  fluence_cm2,
          long*  slab_no,
          float*  parameter,
          long*  N2,
          char*  material_name,
          long*  n_slabs,
          long*  gamma_model,
          float*  gamma_parameter)
{
  long   i;

  FILE*  inputFile;
  inputFile  =  fopen(  fileName, "r");
  char * retval;
  if (inputFile == 0) return;

    char line[256];

    while(fgets(line, 256, inputFile))
    {
    ////////////////////
    // read in parameter
    ////////////////////
    if(strstr(line, ":start parameter:") != NULL){
      while(strstr(line, ":stop parameter:") == NULL)
      {
        retval = fgets(line, 256, inputFile);
        // isolate value from name
        char  seps[]  = "=";                  // delimiters can be ",", ";", or tab
        char*  token;

        token  =  strtok( line, seps );

        if(strstr(token, "r.min.m") != NULL){
          token      =  strtok( NULL, seps );
          parameter[0]  =  (float)atof(token);
        }

        if(strstr(token, "N2") != NULL){
          token    =  strtok( NULL, seps );
          *N2      =  atol(token);
        }


        if(strstr(token, "material.name") != NULL){
          token      =  strtok( NULL, seps );
          char  tmp[256];
          strcpy(tmp, token);

          char  seps2[]  = "\"";                  // delimiter are now '"'
          char*  token2;

          token2      =  strtok( tmp,  seps2);
          token2      =  strtok( NULL, seps2);
          strcpy(material_name, token2);
        }

      }

    }

    ////////////////////
    // read in gamma
    ////////////////////
    if(strstr(line, ":start gamma:") != NULL){
      while(strstr(line, ":stop gamma:") == NULL)
      {
        retval = fgets(line, 256, inputFile);
        // isolate value from name
        char  seps[]  = "=,";                  // delimiters can be ",", ";", or tab
        char*  token;

        token  =  strtok( line, seps );

        if(strstr(token, "model") != NULL){
          token    =  strtok( NULL, seps );
          *gamma_model=  atol(token);

          if(*gamma_model == 1){  *n_gamma_parameter = 2;}
          if(*gamma_model == 2){  *n_gamma_parameter = 3;}
        }

        if(strstr(token, "parameter") != NULL){
          for (i = 0; i < *n_gamma_parameter; i++){
            token        =  strtok( NULL, seps );
            gamma_parameter[i]  =  (float)atof(token);}
        }

      }

    }

    ////////////////////
    // read in spectrum
    ////////////////////
    if(strstr(line, ":start spectrum:") != NULL){
      for (i = 0; i < *nLines; i++){
        retval = fgets(line, 256, inputFile);
        // isolate energy, particle no., fluence, slab no.
        char  seps[]  = ",;\t";                  // delimiters can be ",", ";", or tab
        char*  token;

        token = strtok( line, seps );
        E_MeV_u[i]    =  (float)atof(token);

        token = strtok( NULL, seps );
        particle_no[i]  =  atol(token);

        token = strtok( NULL, seps );
        fluence_cm2[i]  =  (float)atof(token);            // if < 0 --> D.set.Gy

        token = strtok( NULL, seps );
        slab_no[i]    =  atol(token);
      }
    }

    }

  //////////////////////
  // get number of slabs
  //////////////////////

  *n_slabs  =  1;

  for (i = 0; i < *nLines; i++){
    *n_slabs  =  lmaxl(*n_slabs, slab_no[i]);
  }

  return;
}

void AT_readInputS(char**  fileName,
          long*  nLines,
          long*  n_gamma_parameter,
          // return:
          float*  E_MeV_u,
          long*  particle_no,
          float*  fluence_cm2,
          long*  slab_no,
          float*  parameter,
          long*  N2,
          char**  material_name,
          long*  n_slabs,
          long*  gamma_model,
          float*  gamma_parameter){
  AT_readInput(  *fileName,
          nLines,
          n_gamma_parameter,
          // return:
          E_MeV_u,
          particle_no,
          fluence_cm2,
          slab_no,
          parameter,
          N2,
          *material_name,
          n_slabs,
          gamma_model,
          gamma_parameter);};

void AT_browseSpectrum(  char*  fileName,
              // return:
              long*  nLines){
  *nLines = 0;

  FILE*  inputFile;
  inputFile  =  fopen(  fileName, "r");

  if (inputFile == 0) {
    *nLines = -1;
    return;
  }

    char line[256];

  while (fgets(line, 256, inputFile) != NULL){
    *nLines    +=  1;}

  return;
}

void AT_browseSpectrumS(  char**  fileName,
              // return:
              long*  nLines){
  AT_browseSpectrum(  *fileName,
            // return:
            nLines);
};

void AT_readSpectrum(  char*  fileName,
            long*  nLines,
            // return:
            float*  E_MeV_u,
            long*  particle_no,
            float*  fluence_cm2,
            long*  slab_no){
  long   i = 0;

  FILE*  inputFile;
  inputFile  =  fopen(  fileName, "r");

  if (inputFile == 0) return;

    char line[256];

  for (i = 0; i < *nLines; i++){
    if(fgets(line, 256, inputFile) != NULL){
      // isolate energy, particle no., fluence, slab no.
      char  seps[]  = ",;\t";                  // delimiters can be ",", ";", or tab
      char*  token;

      token = strtok( line, seps );
      E_MeV_u[i]    =  (float)atof(token);

      token = strtok( NULL, seps );
      particle_no[i]  =  atol(token);

      token = strtok( NULL, seps );
      fluence_cm2[i]  =  (float)atof(token);            // if < 0 --> D.set.Gy

      token = strtok( NULL, seps );
      slab_no[i]    =  atol(token);
    }
  }
}

void AT_readSpectrumS(  char**  fileName,
            long*  nLines,
            // return:
            float*  E_MeV_u,
            long*  particle_no,
            float*  fluence_cm2,
            long*  slab_no){
  AT_readSpectrum(  *fileName,
            nLines,
            // return:
            E_MeV_u,
            particle_no,
            fluence_cm2,
            slab_no);
};
