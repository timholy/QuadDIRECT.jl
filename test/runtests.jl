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

@testset "Quadratic fit" begin
    a, b, c = rand(), rand(), rand()
    q(x) = a*x^2 + b*x + c
    x1, x2, x3 = sort(rand(3))
    xvert, fvert, qcoef = QuadDIRECT.qfit(x1=>q(x1), x2=>q(x2), x3=>q(x3))
    @test qcoef ≈ a
    @test xvert ≈ -b/(2a)
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
    @test [x.w for x in mes[1]] == [1.0, Inf]
    @test [x.f for x in mes[1]] == [0.0, camel([-2.0,0.0])]
    @test [x.w for x in mes[2]] == [Inf]
    @test [x.f for x in mes[2]] == [0.0]
end

@testset "Sweeps" begin
    splits = ([-2, 0, 2], [-1, 0, 1])
    lower, upper = [-2.75, -1.9], [3.0, 2.0]
    box, x0, xstar = QuadDIRECT.init(camel, splits, lower, upper)
    r1 = QuadDIRECT.get_root(box)
    QuadDIRECT.sweep!(r1, camel, x0, splits, lower, upper)

    splits = ([-2, 0, 2], [-1, 0, 1])
    lower, upper = [-Inf, -Inf], [Inf, Inf]
    box, x0, xstar = QuadDIRECT.init(camel, splits, lower, upper)
    r2 = QuadDIRECT.get_root(box)
    QuadDIRECT.sweep!(r2, camel, x0, splits, lower, upper)

    mn1, mx1 = extrema(r1)
    mn2, mx2 = extrema(r2)
    @test mn1 < mn2
end
