using QuadDIRECT
using Base.Test

@testset "MinimumEdge" begin
    me = MELink{Float64,Float64}('z')
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

    me = MELink{Float64,Float64}(0)
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

