## χ = 16 data, +++-| init
#### path = "./data/chi16"
    g = [0.5, 1.0, 2.0]; beta = [i for i in range(1.,20.,step=0.1)]

## χ = 8 data, ++-| init
#### path = "./data/g_%.1f.jld" g; key(ωn): [f,  ψ _array,  optim]
#### path = "./data/f_and_sx_g_%.1f.txt" g; [β, f_exact, f_chi8, f_chi28, sx_exact, sx_chi8]
    g = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0]
    beta = [i for i in range(1.,40.,step = 0.1)]

#### path = "./data/ug_%.1f.jld" g
    g = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0]
    beta = [10.0, 20.0, 30.0, 40.0]

#### path = "./data/lowT" :low temperature data
    g = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.5, 2.0, 4.0, 6.0, 8.0, 10.0]
    beta = [i for i in range(10,40.,step = 0.1)]

#### path = "./data/gt/gt8_g_%.1f_b_%i.txt" g β : χ = 8
#### path = "./data/gt/gt28_g_%.1f_b_%i.txt" g β : χ = 2*8
    data： ωn  C_mn/(Em - En)（real, imag)
    g = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    number of Masubara frequencies = 20

#### path = "./data/gt/gf_w8_g_%.1f_b_%i.txt" g β : chi = 8
#### path = "./data/gt/gf_w28_g_%.1f_b_%i.txt" g β : chi = 2*8
    data： ωn  G(iωn)（real, imag)
    g = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    number of Masubara frequencies = 20

#### path = "../data/gt/gf_t8_g_%.1f_b_%i.txt" g β : chi = 8
#### path = "../data/gt/gf_t28_g_%.1f_b_%i.txt" g β : chi = 2*8
    data: tau/beta gt
    g = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    beta = [10.0, 20.0, 30.0, 40.0]
    tau = [i for i in range(0, β, length = 100)]



