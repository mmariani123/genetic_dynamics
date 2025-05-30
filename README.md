# ctmc_project

<ins>*Currently under development.</ins> 

The Continuous Time Markov Chain Model (CTMC) Parasite Life Cycle Project - Incrporating the Gillespie Stochastic Simulation Alogirithm (SSA) at certain key steps to generate solutions to the germane differntial equations. Currently uder devleopment (Core code cannot be made pubic currently) and just branched into a second, similar cross-population (both human and vector) disease allele transmission project.

Incorporating MATLAB/OCTAVE, and Rust - for speed critical components (which I have been enjoying learning and deploying) with some Python3 and C code along with parralelization practices and preparation for deoployment on HPCs (High Performance Computer Clusters). I am also using this project to work on my continuous integration knowledge/capabilities. The overal aim of the project is to create a comprehensive program for modelling the plasmodium falcifiparum (malaria parasite) life cycle and track certain alleles of interest down many generations of the mosquito vector. For a similar more basic example germane to this pooject check out resources where people have developed SSA approaches to the SIR (succeptable, infection, recovery) epidemiology model of disease. 

Here is example output from the Gillespie SSA basic SIR modelimplmeneted with parallel processing (100 simulations) in Rust. 

![rust_ssa_sims](https://github.com/user-attachments/assets/29e01851-3221-481e-a5f6-d04761156edb)

The basic example user interface code (Rust, Python3, Shiny for Python) and a few generic functions (above and below) are currently public.

The Previous Shiny for Python example that I created can be found at: https://michael-p-mariani.shinyapps.io/qbs_app/

![genetic_dynamics_combined_user_interfaces_pic](https://github.com/user-attachments/assets/6f5ff517-7101-432a-bc1b-9fa8a7a8f747)

