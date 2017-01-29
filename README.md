#Main files in the respective folders:

CST:  CST.m
4 Noded Elements: Quad.m
8 Noded Elements: Q8_node.m

#The codes require the following:

- Node info (in csv file ): Node ID, X-Coordinate, Y-Coordinate
- Element info (in csv file ): Element ID, Material ID, Element connectivity matrix
- Material info (in txt file ): Material ID, Elastic Modulus, Poisson’s Ratio 

###Boundary Conditions (To be edited in the code):
	– Dirichlet BC: Node ID, DOF, Value
	– Neumann BC: Element ID, Nodes, DOF, Value

Results can be found in the deirectory "Results"

"Report.pdf" contains the report of the project.
