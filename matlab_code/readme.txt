please modify the following two functions to change the shape of the domain. For the time being, the program supports only two kinds of
domain: square & disk, but you can rotate it

get_domain()
is_domain_disk()

To run the program, use the command

fdm_sh_multigrid_2d( 16, 3, -400 )
% base resolution of the grid 16; the multigrid will go to depth 3, i.e., the finest resolution is 64; the lambda = -400

fdm_sh_multigrid_2d( 64, 1, -1000 )
% base resolution of the grid 64; the multigrid will go to depth 1, i.e., the finest resolution is 64; the lambda = -1000
