## chi = 16 data, +++-| init
path = "./data/chi16"
    g = [0.5, 1.0, 2.0]; beta = [i for i in range(1.,20.,step=0.1)]

## chi = 8 data, ++-| init
path = "./data/g_%.1f.jld" g; key(beta): f psi optim
path = "./data/f_and_sx_g_%.1f.txt" g; [beta f_exact f_chi8 f_chi28 sx_exact sx_chi8]
    g = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0]
    beta = [i for i in range(1.,40.,step = 0.1)]

path = "./data/ug_%.1f.jld" g
    g = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0]
    beta = [10.0, 20.0, 30.0, 40.0]

path = "./data/lowT" :low temperature data
    g = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.5, 2.0, 4.0, 6.0, 8.0, 10.0]
    beta = [i for i in range(10,40.,step = 0.1)]

