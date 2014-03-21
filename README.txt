1. What is the effect of the Ks and Kd parameters on the jello? 
If Ks is set to 0, there is no spring between the particles. Ks is the spring's stiffness. By increasing Ks, the spring becomes more stiff and stretches less. By decreasing Ks, the spring becomes less stiff and stretches more.
Kd's effect on the jello is to dissipate the energy of the springs. If Kd is large, the damping force is large, the energy goes down more quickly. If Kd is small, the energy goes down slowly.


2. What are the benefits and the drawbacks of the collision system used here? What are some different ways in which it could be improved?
The collision makes the system more stable in some circumstances. It also simulates the real world in an easier way and produces more realistic scene. For example in the intersection of cylinders, the collision system prevents the jello cube from being pierced through the cylinders. Because the collision system isn't real, it have some drawbacks. Sometimes, it makes the jello easier to explode. And when the particles are inside the layer of collision, there will always be a virtual spring force exerted on it, pushing the particle outside the layer. So if the Ks and Kd parameters are not tuned perfectly, the jello may never stop bouncing because of the collision.
To improve the system, we could add the feature when the jello stops moving, the ks for the virtual spring in the collision becomes 0. We can also implement a more realistic collision system, which might be complicate, but simulates the real world.


3. From lecture, What is an example of a stiff constrained system?
A jello cube, cloth, hair, facial animation.


4. From lecture, What is the difference between and explicit and implicit integration scheme?
For explicit scheme, it calculates the slope by quantities that are known. It requires very small steps for stable integration. For implicit scheme, it calculates slope using y(n+1), which is unknown. And Implicit Euler is stable for all h>0.  


5. Does the jello behave realistically? 
By tuning the parameters perfectly, the jello could bahave realistically. But there're times when it behaves not realistic. For example when a real jello falls on something sharp it will be pierced through and may fall apart. But in the simulation the jello explodes easily.




Here's what I implemented in the project code:

1. Forward Euler and midpoint integration:
   JelloMesh::EulerIntegration()
   JelloMesh::MidpointIntegration()

2. Particle forces other than gravity:
   JelloMesh::ComputeForces(ParticleGrid& grid)

3. Collision and penetration detection:
   JelloMesh::FloorIntersection()
   JelloMesh::CylinderIntersection()
   JelloMesh::SphereIntersection()
   JelloMesh::CubeIntersection()

4. Collision and penetration response:
   JelloMesh::ResolveContacts(ParticleGrid& grid)
   JelloMesh::ResolveCollisions(ParticleGrid& grid)

5. Extra springs:
   I modified JelloMesh::InitJeloMesh() to augment the structural springs with bend and shear springs. 
   There are two types of bend strings. One is to connect every other particles. The second is to connect every two particles.
   I also add diagonal springs that are put on the four diagonals of a particle cube(that is formed by 8 particles), which turned out to make the animation of the jello in the scene of cylinder and sphere to be much more stable.

6.  Support collision with a cube:
   I added a new scene of a cube and implemented JelloMesh::CubeIntersection()

7.  Support collisions with a sphere:
   I also added another new scene of a sphere and implemented JelloMesh::SphereIntersection()
  

8. Tuning the parameters:
   Tuning the parameters turned out to be a tough work. At first I made different sets of parameters for different shapes and different integration methods. After taking plenty of time I finally found three sets of Ks and Kd parameters suitable for RK4, Euler, Midpoint respectively, and one set of the threshold parameters for all shapes and integration methods.


How to run the animation for different integration methods:
For each of the integration method, I have a config file to store the Ks and Kd parameters, which are
RK4: configRK4.txt
EULER: configEULER.txt
MIDPOINT: configMIDPOINT.txt
So if you want to test different methods, you could change the Command Arguments to one of them.
You can also modify the following code in the main.cpp to change between different shapes 
//World theWorld("worlds/ground.xml");
//World theWorld("worlds/cylinders.xml");
//World theWorld("worlds/cube.xml");
//World theWorld("worlds/sphere.xml");


Videos are included in the folder "Output_Videos".

If you have any questions, feel free to contact me: zimengy@seas.upenn.edu








