using QuadDIRECT, StaticArrays
using Base.Test

@testset "MinimumEdge" begin
    me = QuadDIRECT.MELink{Float64,Float64}('z')
    insert!(me, 1, 'a'=>1)
    insert!(me, 2, 'b'=>5)
    @test [x.w for x in me] == [1, 2]
    @test [x.f for x in me] == [1, 5]
    @test [x.l for x in me] == ['a', 'b']
    insert!(me, 4, 'c'=>3)
    @test [x.w for x in me] == [1, 4]
    @test [x.f for x in me] == [1, 3]
    @test [x.l for x in me] == ['a', 'c']
    insert!(me, 10, 'd'=>0)
    @test [x.w for x in me] == [10]
    @test [x.f for x in me] == [0]
    @test [x.l for x in me] == ['d']
    insert!(me, 10, 'e'=>1)
    @test [x.w for x in me] == [10]
    @test [x.f for x in me] == [0]
    @test [x.l for x in me] == ['d']
    insert!(me, 10, 'f'=>0)
    @test [x.w for x in me] == [10]
    @test [x.f for x in me] == [0]
    @test [x.l for x in me] == ['d']
    insert!(me, 10, 'g'=>-1)
    @test [x.w for x in me] == [10]
    @test [x.f for x in me] == [-1]
    @test [x.l for x in me] == ['g']
    insert!(me, 9, 'h'=>-1)
    @test [x.w for x in me] == [10]
    @test [x.f for x in me] == [-1]
    @test [x.l for x in me] == ['g']
    insert!(me, 11, 'i'=>-1)
    @test [x.w for x in me] == [11]
    @test [x.f for x in me] == [-1]
    @test [x.l for x in me] == ['i']

    me = QuadDIRECT.MELink{Float64,Float64}(0)
    y = rand(20)
    w = rand(20)
    for i = 1:20
        insert!(me, w[i], i=>y[i])
    end
    l = sortperm(w)
    w = w[l]
    y = y[l]
    state = start(me)
    i = 1
    ntests = 0
    while !done(me, state)
        item, state = next(me, state)
        while i <= length(w) && w[i] <= item.w
            ntests += 1
            @test y[i] >= item.f
            if y[i] == item.f
                @test l[i] == item.l
            end
            i += 1
        end
    end
    @test ntests > 1
    @test i == 21
end

@testset "Rays" begin
    bb = [(-1,1), (-1,1)]
    t, exitdim = QuadDIRECT.pathlength_box_exit([0, 0.5], [1, 1], bb)
    @test t ≈ 0.5 && exitdim == 2
    t, exitdim = QuadDIRECT.pathlength_box_exit([0, 0.5], [1, 0.3], bb)
    @test t ≈ 1 && exitdim == 1
    t, exitdim = QuadDIRECT.pathlength_box_exit([0, 0.5], [1, -1], bb)
    @test t ≈ 1 && exitdim == 1
    t, exitdim = QuadDIRECT.pathlength_box_exit([0, 0.5], [1, -3], bb)
    @test t ≈ 1/2 && exitdim == 2
    t, exitdim = QuadDIRECT.pathlength_box_exit([0, 0.5], [1, 0], bb)
    @test t ≈ 1 && exitdim == 1
    t, exitdim = QuadDIRECT.pathlength_box_exit([0, 0.5], [0, 0.25], bb)
    @test t ≈ 2 && exitdim == 2
    bb = [(-1,Inf), (-1,1)]
    t, exitdim = QuadDIRECT.pathlength_box_exit([0, 0.5], [1, 0.3], bb)
    @test t ≈ 5/3 && exitdim == 2

    t, intersectdim = QuadDIRECT.pathlength_hyperplane_intersect([0, 0.5], [1, 1], [2, 2], Inf)
    @test t ≈ 2 && intersectdim == 1
    t, intersectdim = QuadDIRECT.pathlength_hyperplane_intersect([0, 0.5], [1, 1], [2, 2], 1.8)
    @test t ≈ 1.5 && intersectdim == 2
    t, intersectdim = QuadDIRECT.pathlength_hyperplane_intersect([0, 0.5], [1, 1], [2, 2], 1)
    @test t == 0 && intersectdim == 0
end

function camel(x)
    # 6-hump camel function. Typically evaluated over [-3,3] × [-2,2].
    x1, x2 = x[1], x[2]
    x1s = x1*x1
    x2s = x2*x2
    return (4 - 2.1*x1s + x1s*x1s/3)*x1s + x1*x2 + (-4 + 4*x2s)*x2s
end

function canyon(x)
    x1, x2 = x[1], x[2]
    return 0.1*(x1+x2)^2 + 10*(x1-x2)^2
end

@testset "Tree topology and parsing" begin
    io = IOBuffer()
    for (str, dim) in (("1(l, l, l)", 1),
                       ("2(l, l, l)", 2))
        box = parse(Box{Float64,3}, str)
        @test QuadDIRECT.isroot(box)
        @test !QuadDIRECT.isleaf(box)
        @test box.splitdim == dim
        for i = 1:3
            @test !QuadDIRECT.isroot(box.children[i])
            @test QuadDIRECT.isleaf(box.children[i])
        end
        splitprint(io, box)
        @test String(take!(io)) == str
    end
    str = "2(l, 1(l, 3(l, l, l), l), l)"
    box = parse(Box{Float64,3}, str)
    splitprint(io, box)
    @test String(take!(io)) == str
    str = "2(l, 1(l, 3(l, l, l), 1(l, l, l)), l)"
    box = parse(Box{Float64,3}, str)
    splitprint(io, box)
    @test String(take!(io)) == str
end

@testset "Initialization" begin
    splits = ([-2, 0, 2], [-1, 0, 1])
    box, x0, xstar = QuadDIRECT.init(camel, splits, [-Inf, -Inf], [Inf, Inf])
    @test x0 == [0, 0]
    @test box isa QuadDIRECT.Box{Float64}
    @test QuadDIRECT.isleaf(box)
    @test !QuadDIRECT.isroot(box)
    p = box.parent
    @test !QuadDIRECT.isleaf(p)
    @test !QuadDIRECT.isroot(p)
    @test p.minmax == (-Inf, Inf)
    @test length(p.children) == 3
    x = [NaN, NaN]
    for i = 1:3
        @test p.xvalues[i] == splits[2][i]
        QuadDIRECT.position!(x, p.children[i])
        @test x == [splits[1][2], splits[2][i]]  # for this initialization, along 1st dim the minimum is in the middle
        @test p.fvalues[i] == camel(x)
    end
    p = p.parent
    @test !QuadDIRECT.isleaf(p)
    @test QuadDIRECT.isroot(p)
    @test p.minmax == (-Inf, Inf)
    @test length(p.children) == 3
    for i = 1:3
        @test p.xvalues[i] == splits[1][i]
        @test p.fvalues[i] == camel([splits[1][i], splits[2][2]]) # perforce the 2nd coordinate is chosen from the middle
    end

    h(x) = sum(abs2, x)/2
    splits = [[-3,-2,-1] for i = 1:3]
    upper = fill(Inf, 3)
    lower = -upper
    box, x0, xstar = QuadDIRECT.init(h, splits, lower, upper)
    root = QuadDIRECT.get_root(box)
    for leaf in leaves(root)
        @test value(leaf) == h(position(leaf, x0))
    end
end

@testset "Traversal" begin
    splits = ([-2, 0, 2], [-1, 0, 1])
    box, x0, xstar = QuadDIRECT.init(camel, splits, [-Inf, -Inf], [Inf, Inf])
    # Finding neighbors
    function nbrprep(child)
        x = QuadDIRECT.position(child, x0)
        i = child.parent.splitdim
        bb = QuadDIRECT.boxbounds(child)
        x, bb, i
    end
    root = QuadDIRECT.get_root(box)
    for box in (root, filter(x->!QuadDIRECT.isleaf(x), root.children)[1])
        x, bb, i = nbrprep(box.children[1])
        x[i] = bb[2]
        nbr, success = QuadDIRECT.find_leaf_at_edge(root, x, i, +1)
        @test QuadDIRECT.isleaf(nbr)
        @test QuadDIRECT.isparent(box.children[2], nbr)
        x, bb, i = nbrprep(box.children[2])
        x[i] = bb[1]
        nbr, success = QuadDIRECT.find_leaf_at_edge(root, x, i, -1)
        @test QuadDIRECT.isleaf(nbr)
        @test QuadDIRECT.isparent(box.children[1], nbr)
        x[i] = bb[2]
        nbr, success = QuadDIRECT.find_leaf_at_edge(root, x, i, +1)
        @test QuadDIRECT.isleaf(nbr)
        @test QuadDIRECT.isparent(box.children[3], nbr)
        x, bb, i = nbrprep(box.children[3])
        x[i] = bb[1]
        nbr, success = QuadDIRECT.find_leaf_at_edge(root, x, i, -1)
        @test QuadDIRECT.isleaf(nbr)
        @test QuadDIRECT.isparent(box.children[2], nbr)
    end
end

@testset "Quadratic fit" begin
    # qfit tests
    a, b, c = rand(), rand(), rand()
    q(x) = a*x^2 + b*x + c
    x1, x2, x3 = sort(rand(3))
    xvert, fvert, qcoef = QuadDIRECT.qfit(x1=>q(x1), x2=>q(x2), x3=>q(x3))
    @test qcoef ≈ a
    @test xvert ≈ -b/(2a)

    # qdelta tests
    q(x, y, z) = x^2 - y^2 + z
    x0 = y0 = z0 = 0
    root = QuadDIRECT.Box{Float64,3}()
    QuadDIRECT.add_children!(root, 1, [-1, 0.2, 1], [q(-1,y0,z0), q(0.2,y0,z0), q(1,y0,z0)], -Inf, Inf)
    @test QuadDIRECT.qdelta(root.children[1]) ≈ q(-0.4,y0,z0) - q(-1,y0,z0)
    @test QuadDIRECT.qdelta(root.children[2]) ≈ q(0,y0,z0) - q(0.2,y0,z0)
    @test QuadDIRECT.qdelta(root.children[3]) ≈ q(0.6,y0,z0) - q(1,y0,z0)
    for i = 1:3
        @test QuadDIRECT.qdelta(root.children[i], 1) == QuadDIRECT.qdelta(root.children[i])
        @test QuadDIRECT.qdelta(root.children[i], 2) == 0
        @test QuadDIRECT.qdelta(root.children[i], 3) == 0
    end
    p = root.children[3]
    x0 = p.parent.xvalues[p.parent_cindex]
    QuadDIRECT.add_children!(p, 2, [-1, 0.2, 1], [q(x0,-1,z0), q(x0,0.2,z0), q(x0,1,z0)], -Inf, Inf)
    @test QuadDIRECT.qdelta(p.children[1]) == -Inf
    @test QuadDIRECT.qdelta(p.children[2]) ≈ q(x0,0.6,z0) - q(x0,0.2,z0)
    @test QuadDIRECT.qdelta(p.children[3]) == -Inf
    for i = 1:3
        @test QuadDIRECT.qdelta(p.children[i], 2) == QuadDIRECT.qdelta(p.children[i])
        @test QuadDIRECT.qdelta(p.children[i], 1) == QuadDIRECT.qdelta(p)
        @test QuadDIRECT.qdelta(p.children[i], 3) == 0
    end
    p = p.children[1]
    y0 = p.parent.xvalues[p.parent_cindex]
    QuadDIRECT.add_children!(p, 3, [-1, 0.2, 1], [q(x0,y0,-1), q(x0,y0,0.2), q(x0,y0,1)], -Inf, Inf)
    @test QuadDIRECT.qdelta(p.children[1]) == -Inf
    @test QuadDIRECT.qdelta(p.children[2]) ≈ q(x0,y0,-0.4) - q(x0,y0,0.2)
    @test QuadDIRECT.qdelta(p.children[3]) ≈ q(x0,y0, 0.6) - q(x0,y0,1)
    for i = 1:3
        @test QuadDIRECT.qdelta(p.children[i], 3) == QuadDIRECT.qdelta(p.children[i])
        @test QuadDIRECT.qdelta(p.children[i], 2) == QuadDIRECT.qdelta(p)
        @test QuadDIRECT.qdelta(p.children[i], 1) == QuadDIRECT.qdelta(p.parent)
    end

    # quasinewton tests
    splits = ([-11,-10,-9], [-7,-6,-5])
    lower, upper = [-Inf, -Inf], [Inf, Inf]
    box, x0, xstar = QuadDIRECT.init(canyon, splits, lower, upper)
    Q, xbase, c = QuadDIRECT.build_quadratic_model(box, x0)
    @test xbase == xstar
    @test c == canyon(xstar)
    @test_throws Base.LinAlg.SingularException QuadDIRECT.solve(Q)
    y = xstar[2]
    QuadDIRECT.add_children!(box, 1, [xstar[1], -8, -7],
                             [canyon([xstar[1],y]), canyon([-8,y]), canyon([-7,y])], -Inf, Inf)
    root = QuadDIRECT.get_root(box)
    for leaf in leaves(root)
        Q, xbase, c = QuadDIRECT.build_quadratic_model(leaf, x0)
        g, B = QuadDIRECT.solve(Q)
        @test B ≈ [20.2 -19.8; -19.8 20.2]
        @test (B \ g) ≈ xbase
    end
end

@testset "Building minimum edges" begin
    splits = ([-2, 0, 2], [-1, 0, 1])
    lower, upper = [-2.75, -1.9], [3.0, 2.0]
    box, x0, xstar = QuadDIRECT.init(camel, splits, lower, upper)
    r = QuadDIRECT.get_root(box)
    mes = QuadDIRECT.minimum_edges(r, x0, lower, upper)
    @test [x.w for x in mes[1]] == [1.0]
    @test [x.f for x in mes[1]] == [0.0]
    @test [x.w for x in mes[2]] == [1.0, 2.0]
    @test [x.f for x in mes[2]] == [0.0, camel([-2.0,0.0])]
    splits = ([-2, 0, 2], [-1, 0, 1])
    lower, upper = [-Inf, -Inf], [Inf, Inf]
    box, x0, xstar = QuadDIRECT.init(camel, splits, lower, upper)
    r = QuadDIRECT.get_root(box)
    mes = QuadDIRECT.minimum_edges(r, x0, lower, upper)
    @test [x.w for x in mes[1]] == [1.0]
    fs = [x.f for x in mes[1]]
    bxs = [x.l for x in mes[1]]
    @test length(fs) == 1 && all([fs[i] <= camel(position(bxs[i], x0)) for i = 1:1])
    @test [x.w for x in mes[2]] == [0.5, Inf]
    fs = [x.f for x in mes[2]]
    bxs = [x.l for x in mes[2]]
    @test length(fs) == 2 && all([fs[i] <= camel(position(bxs[i], x0)) for i = 1:2])
end

@testset "Sweeps" begin
    splits = ([-2, 0, 2], [-1, 0, 1])
    lower, upper = [-2.75, -1.9], [3.0, 2.0]
    box, x0, xstar = QuadDIRECT.init(camel, splits, lower, upper)
    r1 = QuadDIRECT.get_root(box)
    QuadDIRECT.sweep!(r1, QuadDIRECT.CountedFunction(camel), x0, splits, lower, upper)

    splits = ([-2, 0, 2], [-1, 0, 1])
    lower, upper = [-Inf, -Inf], [Inf, Inf]
    box, x0, xstar = QuadDIRECT.init(camel, splits, lower, upper)
    r2 = QuadDIRECT.get_root(box)
    QuadDIRECT.sweep!(r2, QuadDIRECT.CountedFunction(camel), x0, splits, lower, upper)

    mn1, mx1 = extrema(r1)
    mn2, mx2 = extrema(r2)
    @test mn1 < mn2
end

@testset "Minimize" begin
    splits = ([-2, 0, 2], [-1, 0, 1])
    lower, upper = [-3,-2], [3,2]
    xmin, fmin = minimize(camel, splits, lower, upper)
    @test fmin < -1.02
    @test norm(xmin - [0.0898, -0.7126]) < 0.1 || norm(xmin - [-0.0898, 0.7126]) < 0.1
    root, x0 = analyze(camel, splits, lower, upper)
    @test length(leaves(root)) < 500

    splits = ([-11,-10,-9], [-7,-6,-5])
    lower, upper = [-Inf, -Inf], [Inf, Inf]
    xmin, fmin = minimize(canyon, splits, lower, upper; atol=1e-3)
    @test fmin < 1e-10
    @test norm(xmin) < 1e-5
    root, x0 = analyze(canyon, splits, lower, upper)
    @test length(leaves(root)) < 300
end

@testset "Infinite return values" begin
    canyonb(x) = x[2] > -0.1 ? Inf : canyon(x)
    splits = ([-11,-10,-9], [-7,-6,-5])
    lower, upper = [-Inf, -Inf], [Inf, Inf]
    root, x0 = analyze(canyonb, splits, lower, upper)
    box = minimum(root)
    x = position(box, x0)
    @test x[2] <= -0.1
    @test value(box) < 0.01
    @test length(leaves(root)) < 500
end

@testset "High dimensional" begin
    B = randn(41, 40); B = B'*B
    D, V = eig(B)
    i = 1
    while D[i] < 1e-4*D[end]
        D[i] = 1e-4*D[end]
    end
    B = V*Diagonal(sqrt.(D)); B = B'*B
    h(x) = (x'*B*x)/2
    splits = [[-3,-2,-1] for i = 1:size(B,1)]
    upper = fill(Inf, size(B,1))
    lower = -upper
    root, x0 = analyze(h, splits, lower, upper; maxevals=10000, fvalue=1e-3)
    @test value(minimum(root)) <= 1e-3
end
