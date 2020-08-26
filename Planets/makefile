# Makefile for Planets N-Body interaction solver

GEOMETRY_DIR = Geometry/
BODY_DIR = Body/
SPH_DIR = SPH/
EULER_DIR = Euler/
MPI_DIR = Planets_MPI/
OCTREE_DIR = Octree/
BINARYTREE_DIR = BinaryTree/
HASHGRID_DIR = SPH/
INPUT_DIR = Input/
MAIN_DIR = Main/

Planets_NBody: Planets_NBody.o Body.o Multi_Body_System.o MathVector.o Planets_MPI.o Octree.o BinaryTree.o Input.o
	mpiCC -o Planets_NBody Planets_NBody.o Body.o Multi_Body_System.o MathVector.o Planets_MPI.o Octree.o BinaryTree.o Input.o

Planets_SPH: Planets_SPH.o BodySPH.o MathVector.o Planets_MPI.o BinaryTree.o HashGrid.o Input.o
	mpiCC -o Planets_SPH Planets_SPH.o BodySPH.o MathVector.o Planets_MPI.o BinaryTree.o Input.o

Planets_NBody.o: $(MAIN_DIR)Planets_NBody.cc
	mpiCC -c -Wall $(MAIN_DIR)Planets_NBody.cc

Planets_SPH.o: $(MAIN_DIR)Planets_SPH.cc
	mpiCC -c -Wall $(MAIN_DIR)Planets_SPH.cc

Planets_Euler.o: $(MAIN_DIR)Planets_Euler.cc
	mpiCC -c -Wall $(MAIN_DIR)Planets_Euler.cc

Body.o: $(BODY_DIR)Body.cc
	icpc -c -Wall $(BODY_DIR)Body.cc

BodySPH.o: $(SPH_DIR)BodySPH.cc
	icpc -c -Wall $(SPH_DIR)BodySPH.cc

Multi_Body_System.o: $(BODY_DIR)Multi_Body_System.cc
	mpiCC -c -Wall $(BODY_DIR)Multi_Body_System.cc

MathVector.o: $(GEOMETRY_DIR)MathVector.cc
	icpc -c -Wall $(GEOMETRY_DIR)MathVector.cc

Planets_MPI.o: $(MPI_DIR)Planets_MPI.cc
	mpiCC -c -Wall $(MPI_DIR)Planets_MPI.cc

Octree.o: $(OCTREE_DIR)Octree.cc
	icpc -c -Wall $(OCTREE_DIR)Octree.cc

BinaryTree.o: $(BINARYTREE_DIR)BinaryTree.cc
	icpc -c -Wall $(BINARYTREE_DIR)BinaryTree.cc

HashGrid.o: $(HASHGRID_DIR)HashGrid.cc
	icpc -c -Wall $(HASHGRID_DIR)HashGrid.cc

Input.o: $(INPUT_DIR)Input.cc
	mpiCC -c -Wall $(INPUT_DIR)Input.cc

clean:
	rm -f *.o
	rm -f *~
	rm -f $(GEOMETRY_DIR)*.o
	rm -f $(GEOMETRY_DIR)*~
	rm -f $(BODY_DIR)*.o
	rm -f $(BODY_DIR)*~
	rm -f $(MPI_DIR)*.o
	rm -f $(MPI_DIR)*~
	rm -f $(OCTREE_DIR)*.o
	rm -f $(OCTREE_DIR)*~
	rm -f $(BINARYTREE_DIR)*.o
	rm -f $(BINARYTREE_DIR)*~
	rm -f $(EULER_DIR)*.o
	rm -f $(EULER_DIR)*~