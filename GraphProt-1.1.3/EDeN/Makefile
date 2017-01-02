# SVM with stochastic gradient

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA

#define HINGELOSS 1
#define SMOOTHHINGELOSS 2
#define SQUAREDHINGELOSS 3
#define LOGLOSS 10
#define LOGLOSSMARGIN 11

#CHECKLIMITS activates the limit check for graph operations

CXX=g++
C=gcc
OPTS=-g -O3 -Wno-deprecated -DLOSS=1  -fopenmp -static  #-DDEBUGON -DCHECKLIMITS
#OPTS=-g -Wno-deprecated -DLOSS=1 -DDEBUGON
#OPTS=-g -Wno-deprecated -static -pg
CXXFLAGS=${OPTS} -Wall -mtune=native -march=native
#OBABEL=-I /usr/local/openbabel/2.3.1/include/openbabel-2.0  -L /usr/local/openbabel/2.3.1/lib -DUSEOBABEL
#OBABEL=-I /usr/include/openbabel-2.0 -DUSEOBABEL
#CXXFLAGS= ${OPTS} 
LIBS=-lm -lz #-lopenbabel 

PROGRAMS=EDeN 

all: ${PROGRAMS}

clean:
	-rm ${PROGRAMS}
	-rm *.o

KQuickShift:

NSPDK: NSPDK.o vectors.o gzstream.o GraphClass.o NSPDK_FeatureGenerator.o vectors.o FlagsService.o Utility.o OpenBabelConverter.o                                            
	${CXX} ${CXXFLAGS} ${OBABEL}  NSPDK.o vectors.o gzstream.o BaseGraphClass.o GraphClass.o NSPDK_FeatureGenerator.o FlagsService.o Utility.o OpenBabelConverter.o ${LIBS} -o $@   

EDeN: EDeN.o vectors.o gzstream.o GraphClass.o NSPDK_FeatureGenerator.o vectors.o FlagsService.o Utility.o OpenBabelConverter.o DDK.o
	${CXX} ${CXXFLAGS} ${OBABEL} EDeN.o vectors.o gzstream.o BaseGraphClass.o GraphClass.o DDK.o NSPDK_FeatureGenerator.o FlagsService.o Utility.o OpenBabelConverter.o ${LIBS} -o $@

EDeN.o: EDeN.cpp vectors.h  gzstream.h OpenBabelConverter.o
	 ${CXX} ${CXXFLAGS} -c EDeN.cpp ${OBABEL} -o EDeN.o

FlagsService.o:FlagsService.h

NSPDK_FeatureGenerator.o:NSPDK_FeatureGenerator.h GraphClass.o 

GraphClass.o:GraphClass.h BaseGraphClass.o 

BaseGraphClass.o:BaseGraphClass.h 

OpenBabelConverter.o:OpenBabelConverter.h GraphClass.h BaseGraphClass.o
	${CXX} ${CXXFLAGS} -c OpenBabelConverter.cpp ${OBABEL} -o OpenBabelConverter.o

Utility.o:Utility.h vectors.h

vectors.o: vectors.cpp vectors.h  wrapper.h

gzstream.o: gzstream.cpp gzstream.h 

DDK.o: DDK_FeatureGenerator.h DDK_FeatureGenerator.cpp NSPDK_FeatureGenerator.o
	$(CXX) $(CFLAGS) -c DDK_FeatureGenerator.cpp $(LIBS) $(OBABEL) -o DDK.o
