import os
import logging
import fipy
from fipy import CellVariable, Grid3D, Grid1D, DiffusionTerm, ConvectionTerm, TransientTerm, CylindricalGrid2D, CentralDifferenceConvectionTerm

# import meshio
log = logging.getLogger("fipy")

console = logging.StreamHandler()

console.setLevel(logging.INFO)

log.addHandler(console)


if __name__ =="__main__":

    # filename = "hollow_cylinder.msh"
    # vertices = []
    # faces = []

    # with open(filename, "r") as f:
    #     reading_nodes = False
    #     reading_elements = False

    #     for line in f:
    #         if line.startswith("$Nodes"):
    #             reading_nodes = True
    #         elif line.startswith("$EndNodes"):
    #             reading_nodes = False
    #         elif line.startswith("$Elements"):
    #             reading_elements = True
    #         elif line.startswith("$EndElements"):
    #             reading_elements = False
    #         elif reading_nodes:
    #             vertex_data = line.strip().split()
    #             if len(vertex_data) == 4:
    #                 vertices.append([float(vertex_data[1]), float(vertex_data[2]), float(vertex_data[3])])
    #         elif reading_elements:
    #             element_data = line.strip().split()
    #             if len(element_data) >= 8 and element_data[1] == "4":  # Tetrahedral elements
    #                 element_vertices = [int(v) - 1 for v in element_data[8:12]]
    #                 faces.append(element_vertices)

    # # Create the FiPy mesh
    # mesh = Grid3D(vertices=vertices, faces=faces)


# Define the dimensions of the hollow cylinder
    inner_radius = 1.0
    outer_radius = 2.0
    length = 5.0
    num_cells = 50

    # Create a cylindrical 2D mesh with the hollow cylinder geometry
    mesh = CylindricalGrid2D(dx=length/num_cells, dy=(outer_radius-inner_radius)/num_cells, nx=num_cells, ny=num_cells)

  # Import the Gmsh MSH mesh into FiPy
    # gmsh_mesh = meshio.read("your_gmsh_mesh.msh")
#   mesh = Grid3D.fromGmsh("hollow_cylinder.msh")

  # Define your FiPy simulation using the imported mesh


    rho = CellVariable(mesh=mesh, name="density", value=0.0)
    u = CellVariable(mesh=mesh, name="velocity", value=(0.0, 0.0))
    p = CellVariable(mesh=mesh, name="pressure", value=1.0)
    tg = CellVariable(mesh=mesh, name="pressure", value=1.0)
    e = CellVariable(mesh=mesh, name="energy", value=0.0)

  # ... Add your FiPy equations, boundary conditions, and solver here ...


# Define all equation sets using terms

    Q = 100.0  # Heat source (can be zero if there is no heat source)
    Cp = 1000.0  # Specific heat capacity at constant pressure
    k = 0.1  # Thermal conductivity

# Define the diffusion coefficient
    D = 1.0
# Define the convection velocity
    v = (1.0,)


# EQUATION FIELD RECONSTRUCTION
# Define the equation with second-order reconstruction using CentralDifferenceConvectionTerm
    eq_continuity = TransientTerm() == DiffusionTerm(coeff=1.0) - CentralDifferenceConvectionTerm(coeff=1.0, var=rho) 

# Define the continuity equation
    # eq_continuity = TransientTerm(coeff=rho) == -ConvectionTerm(coeff=u.divergence)

# Define the momentum equations
    eq_momentum_x = (TransientTerm(coeff=rho, var=u[0]) + ConvectionTerm(coeff=rho*u[0], var=u[0])
                == -DiffusionTerm(coeff=p) + rho.faceGrad.divergence)

    eq_momentum_r = (TransientTerm(coeff=rho, var=u[1]) + ConvectionTerm(coeff=rho*u[1], var=u[1])
                == -DiffusionTerm(coeff=p) + rho.faceGrad.divergence)

# Define the energy equation
    eq_energy = (TransientTerm(coeff=rho*Cp, var=tg) + ConvectionTerm(coeff=rho*Cp*u, var=tg)
            == DiffusionTerm(coeff=k, var=T) + Q)


# Apply the limiter to the convection term
    convection_limiter = vanLeerLimiter
    eq._terms[1]._convective._diffusion.limiter = convection_limiter



# Boundary conditions

# Dirichlet
    rho_value_left = 0.0
    rho.constrain(rho_value_left, mesh.facesLeft)


# Neumann
    gradient_right = 1.0
    rho.faceGrad.constrain(gradient_right, mesh.facesRight)


# Apply no-slip boundary condition along the inner surface of the hollow cylinder
    u[0].constrain(0.0, where=mesh.x < inner_radius)
    u[1].constrain(0.0, where=mesh.x < inner_radius)

# Flux Boundary conditions
    flux_value_at_boundary = 2.0
    var.constrainFlux(flux_value_at_boundary, mesh.facesRight)  # Set the flux at the right boundary

# Boundary conditions Velocity
    value_left = 0.0
    rho.constrain(value_left, mesh.facesLeft)

    gradient_right = 1.0
    rho.faceGrad.constrain(gradient_right, mesh.facesRight)

# Boundary conditions Temperatures
    value_left = 0.0
    rho.constrain(value_left, mesh.facesLeft)

    gradient_right = 1.0
    rho.faceGrad.constrain(gradient_right, mesh.facesRight)


# Boundary conditions Outlet free flow - Non-reflective
    value_left = 0.0
    rho.constrain(value_left, mesh.facesLeft)

    gradient_right = 1.0
    rho.faceGrad.constrain(gradient_right, mesh.facesRight)


# Solve the equation
    time_step = 0.001
    total_time = 1.0
    steps = int(total_time / time_step)
    for step in range(steps):
        rho.updateOld()  # Store the current solution as the old solution
        u.updateOld()  # Store the current solution as the old solution
        v.updateOld()  # Store the current solution as the old solution
        e.updateOld()  # Store the current solution as the old solution
        tg.updateOld()  # Store the current solution as the old solution
        p.updateOld()  # Store the current solution as the old solution
        
        eq_continuity.solve(var=rho, dt=time_step)  # Solve the x-component of momentum equation        # Optionally, update the Viewer at each time step
        eq_momentum_x.solve(var=u[0], dt=time_step)  # Solve the x-component of momentum equation        # Optionally, update the Viewer at each time step
        eq_momentum_r.solve(var=u[1], dt=time_step)  # Solve the x-component of momentum equation        # Optionally, update the Viewer at each time step
        eq_energy.solve(var=tg, dt=time_step)  # Solve the energy equation for temperature


# Apply the limiter after each time step to control the solution near discontinuities
        rho.constrain(rho, mesh.exteriorFaces)
        u[0].constrain(u[0], mesh.exteriorFaces)
        u[1].constrain(u[1], mesh.exteriorFaces)
        tg.constrain(tg, mesh.exteriorFaces)


        if step % 10 == 0:
            viewer.plot()


  # Optionally, you can use FiPy's Viewer to visualize the simulation


    viewer = Viewer(vars=phi)
    viewer.plot()

  # Run your FiPy simulation
  # ... (Time-stepping and solving the equations) ...

