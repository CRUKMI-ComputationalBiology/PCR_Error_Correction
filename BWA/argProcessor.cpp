/*
  Copyright (C) 2017, Cancer Research UK Manchester Institute

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "argProcessor.hpp"
#include <stdlib.h>

ArgProcessor::ArgProcessor(int argc, char* argv[]) {
  
  for (int i=1; i < argc; i++) {
	char* arg = argv[i];
	string myArg (arg);
	if (arg[0] == '-') { 
	
	  argSet[myArg] = true;
  	  	  
	  if (i != argc-1) {
		string nextArg (argv[i+1]);
		argVal[myArg] = nextArg;
	  }
	}
  }
  
}


bool ArgProcessor::isArgSet(string arg) {
  
  map<string,bool>::const_iterator it;

  it = argSet.find(arg);
  if (it != argSet.end()) {
	return(true);
  }
  else {
	return(false);
  }
}

int ArgProcessor::getIntVal(string arg) {
  
  return(atoi(argVal[arg].c_str()));
  
}

long ArgProcessor::getLongVal(string arg) {
  return(atol(argVal[arg].c_str()));
}


float ArgProcessor::getFloatVal(string arg) {
  
  return(atof(argVal[arg].c_str()));
}


string ArgProcessor::getStringVal(string arg) {
  
  return(argVal[arg]);
}


