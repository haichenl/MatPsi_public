/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <string>
#include <map>
#include <sstream>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libparallel/parallel.h>

#include <boost/filesystem.hpp>

namespace psi {

void PSIO::open(unsigned int unit, int status) {
    
    // check if our working folder exists 
    if ( !boost::filesystem::exists(_psio_manager_->get_default_path()) ) {
        boost::filesystem::create_directories(_psio_manager_->get_default_path());
    }
    
  unsigned int i;
  char *name, *path;
  psio_ud *this_unit;
  
  //std::cout << "proc " << WorldComm->me() << "    status = " << status << std::endl;
  /* check for too large unit */
  if (unit > PSIO_MAXUNIT)
    psio_error(unit, PSIO_ERROR_MAXUNIT);
  
  this_unit = &(psio_unit[unit]);
  
  /* Get number of volumes to stripe across */
  this_unit->numvols = get_numvols(unit);
  if (this_unit->numvols > PSIO_MAXVOL)
    psio_error(unit, PSIO_ERROR_MAXVOL);
  if (!(this_unit->numvols))
    this_unit->numvols = 1;
  
  /* Check to see if this unit is already open */
  for (i=0; i < this_unit->numvols; i++) {
    if (this_unit->vol[i].stream != -1)
      psio_error(unit, PSIO_ERROR_REOPEN);
  }
  
  /* Get the file name prefix */
  get_filename(unit, &name);
  //printf("%s\n",name);
  
  // Check if any files will have the same name
  {
    using std::string;
    typedef std::map<string,int> Names;
    Names names;
    for (i=0; i < this_unit->numvols; i++) {
      std::ostringstream oss;
      get_volpath(unit, i, &path);
      oss << path << name << "." << unit;
      const std::string fullpath = oss.str();
      typedef Names::const_iterator citer;
      citer n = names.find(fullpath);
      if (n != names.end())
        psio_error(unit, PSIO_ERROR_IDENTVOLPATH);
      names[fullpath] = 1;
      free(path);
    }
  }
  
  /* Build the name for each volume and open the file */
  for (i=0; i < this_unit->numvols; i++) {
    char* fullpath;
    get_volpath(unit, i, &path);

    #pragma warn A bit of a hack in psio open at the moment, breaks volumes and some error checking
    const char* path2 = _psio_manager_->get_file_path(unit).c_str(); 
    
    fullpath = (char*) malloc( (strlen(path2)+strlen(name)+80)*sizeof(char));
    sprintf(fullpath, "%s%s.%u", path2, name, unit);
    this_unit->vol[i].path = strdup(fullpath);
    free(fullpath);
    
    /* Register the file */
    _psio_manager_->open_file(std::string(this_unit->vol[i].path), unit);

    /* Now open the volume */
    if (status == PSIO_OPEN_OLD) {
      //~ if (WorldComm->me() == 0) {
        this_unit->vol[i].stream = ::open(this_unit->vol[i].path,O_CREAT|O_RDWR,0644);
      //~ }
      //~ WorldComm->bcast(&(this_unit->vol[i].stream), 1, 0);
      //WorldComm->raw_bcast(&(this_unit->vol[i].stream), sizeof(int), 0);
      if(this_unit->vol[i].stream == -1)
        psio_error(unit,PSIO_ERROR_OPEN);
    }
    else if(status == PSIO_OPEN_NEW) {
      //~ if (WorldComm->me() == 0) {
        this_unit->vol[i].stream = ::open(this_unit->vol[i].path,O_CREAT|O_RDWR|O_TRUNC,0644);
      //~ }
      //~ WorldComm->bcast(&(this_unit->vol[i].stream), 1, 0);
      //WorldComm->raw_bcast(&(this_unit->vol[i].stream), sizeof(int), 0);
      if(this_unit->vol[i].stream == -1)
        psio_error(unit,PSIO_ERROR_OPEN);
    }
    else psio_error(unit,PSIO_ERROR_OSTAT);

    free(path);
  }

  if (status == PSIO_OPEN_OLD) tocread(unit);
  else if (status == PSIO_OPEN_NEW) {
    /* Init the TOC stats and write them to disk */
    this_unit->toclen = 0;
    this_unit->toc = NULL;
    wt_toclen(unit, 0);
  }
  else psio_error(unit,PSIO_ERROR_OSTAT);

  free(name);
}

// Mirrors PSIO::open() but just check to see if the file is there
// status needs is assumed PSIO_OPEN_OLD if this is called
bool PSIO::exists(unsigned int unit) {
  unsigned int i;
  char *name, *path;
  psio_ud *this_unit;
  
  if (unit > PSIO_MAXUNIT)
    psio_error(unit, PSIO_ERROR_MAXUNIT);
  
  this_unit = &(psio_unit[unit]);
  
  /* Get number of volumes to stripe across */
  this_unit->numvols = get_numvols(unit);
  if (this_unit->numvols > PSIO_MAXVOL)
    psio_error(unit, PSIO_ERROR_MAXVOL);
  if (!(this_unit->numvols))
    this_unit->numvols = 1;
  
  /* Check to see if this unit is already open, if so, should be good.
     If every volume has a sream value other than -1, it's open */
  bool already_open = true;
  for (i=0; i < this_unit->numvols; i++) {
    if (this_unit->vol[i].stream == -1)
      already_open = false;
  }
  if (already_open) return(true);
  
  /* Get the file name prefix */
  get_filename(unit, &name);
  //printf("%s\n",name);
  
  // Check if any files will have the same name
  {
    using std::string;
    typedef std::map<string,int> Names;
    Names names;
    for (i=0; i < this_unit->numvols; i++) {
      std::ostringstream oss;
      get_volpath(unit, i, &path);
      oss << path << name << "." << unit;
      const std::string fullpath = oss.str();
      typedef Names::const_iterator citer;
      citer n = names.find(fullpath);
      if (n != names.end())
        psio_error(unit, PSIO_ERROR_IDENTVOLPATH);
      names[fullpath] = 1;
      free(path);
    }
  }
  
  /* Build the name for each volume and open the file */
  bool file_exists = true;
  for (i=0; i < this_unit->numvols; i++) {
    char* fullpath;
    int stream;
    get_volpath(unit, i, &path);

    #pragma warn A bit of a hack in psio open at the moment, breaks volumes and some error checking
    const char* path2 = _psio_manager_->get_file_path(unit).c_str(); 
    
    fullpath = (char*) malloc( (strlen(path2)+strlen(name)+80)*sizeof(char));
    sprintf(fullpath, "%s%s.%u", path2, name, unit);
    
    /* Now open the volume */
    //~ if (WorldComm->me() == 0) {
      stream = ::open(fullpath,O_RDWR);
    //~ }
    if (stream == -1) {
      file_exists = false;
    }

    free(path);
    free(fullpath);
  }

  free(name);
  return(file_exists);
}



void
PSIO::rehash(unsigned int unit)
{
  if (open_check(unit)) {
    close(unit,1);
    open(unit,PSIO_OPEN_OLD);
  }
}

  int psio_open(unsigned int unit, int status) {
    _default_psio_lib_->open(unit, status);
    return 1;
  }   

}

