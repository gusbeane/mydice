[1]: https://bitbucket.org/vperret/dice/wiki/browse/
[2]: https://bitbucket.org/rteyssie/ramses
[3]: https://bitbucket.org/rteyssie/ramses/src/fc657574bb4a06b3647d178a88cfe7943b02a4fd/trunk/ramses/patch/init/dice/?at=master

## This is the [DICE](https://bitbucket.org/vperret/dice) bitbucket repository.

DICE is an open source code modelling initial conditions of idealised galaxies to study their secular evolution, or to study more complex interactions such as mergers or compact groups using N-Body/hydro codes.
The particularity of this code is its ability to setup a large number of components modelling distinct parts of the galaxy.
The code creates 3D distributions of particles using a N-try MCMC algorithm which does not require a prior knowledge of the distribution function. The gravitational potential is then computed on a multi-level cartesian mesh by solving the poisson equation in the Fourier space.
Finally, the dynamical equilibrium of each component is computed by integrating the Jeans equations for each particles.
Several galaxies can be generated in a row and be placed on Keplerian orbits to model interactions.
DICE writes the initial conditions in the **Gadget1** or **Gadget2** format and is fully compatible with [Ramses][2] thanks to a [patch][3] included in the public ramses distribution.

The user's guide is accessible [here][1].

Download the code by cloning the git repository using 
```
$ git clone https://bitbucket.org/vperret/dice
```
Please register also to the [DICE google group](http://groups.google.com/forum/#!forum/dice-project).