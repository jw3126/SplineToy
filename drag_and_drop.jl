using GLMakie
import MakieCore as MC
using ArgCheck
using StaticArrays
using OffsetArrays

function make_drag_and_drop()
    positions = Observable([Point2f(x, x) for x in range(0, 1, length=2)])
    dragging = false
    idx = 1

    fig, ax, p = scatter(positions, color=:gold)

    on(events(fig).mousebutton, priority = 2) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                plt, i = pick(fig)
                if Keyboard.d in events(fig).keyboardstate && plt == p
                    # Delete marker
                    deleteat!(positions[], i)
                    notify(positions)
                    return Consume(true)
                elseif Keyboard.a in events(fig).keyboardstate
                    # Add marker
                    push!(positions[], mouseposition(ax))
                    notify(positions)
                    return Consume(true)
                else
                    # Initiate drag
                    dragging = plt == p
                    idx = i
                    return Consume(dragging)
                end
            elseif event.action == Mouse.release
                # Exit drag
                dragging = false
                return Consume(false)
            end
        end
        return Consume(false)
    end

    on(events(fig).mouseposition, priority = 2) do mp
        if dragging
            positions[][idx] = mouseposition(ax)
            notify(positions)
            return Consume(true)
        end
        return Consume(false)
    end

    (;fig, ax, positions)
end

################################################################################
#### Hermite
################################################################################

struct Hermite
    points::Vector{SVector{2,Float32}}
end

function tlims(h::Hermite)
    (0, length(h.points)-1)
end

function estiamte_tangent(points, i)::eltype(points)
    if i == firstindex(points)
        return points[i+1] - points[i]
    elseif i == lastindex(points)
        return points[i] - points[i-1]
    else
        return (points[i+1] - points[i-1]) / 2
    end
end

function (h::Hermite)(t)::SVector{2,Float32}
    t_min, t_max = tlims(h)
    @argcheck t_min <= t <= t_max
    if t == t_max
        return h.points[end]
    end
    n = length(h.points)
    istart = floor(Int, t) + firstindex(h.points)
    istop = istart+1
    
    t = t - floor(t)
    @check 0 <= t <= 1
    h00 = (1 + 2*t) * (1-t)^2
    h10 = t * (1-t)^2
    h01 = t^2 * (3-2*t)
    h11 = t^2 * (t-1)

    pt_start = h.points[istart]
    pt_stop = h.points[istop]
    tangent_start = estiamte_tangent(h.points, istart)
    tangent_stop = estiamte_tangent(h.points, istop)
    
    h00 * pt_start + h10 * tangent_start + h01 * pt_stop + h11 * tangent_stop
end

function MC.plottype(::Hermite)
    MC.Lines
end

function MC.convert_arguments(trait::MC.PointBased, h::Hermite)
    t_min, t_max = tlims(h)
    pts = Point2f[h(t) for t in range(t_min, t_max, length=1000)]
    MC.convert_arguments(trait, pts)
end

################################################################################
#### BSpline
################################################################################
struct BSpline
    knots::Vector{Float64}
    index::Int
    degree::Int
    function BSpline(knots, index, degree)
        @argcheck degree >= 0
        @argcheck issorted(knots)
        @argcheck index in 1:nbsplines(;nknots=length(knots), degree)
        new(knots, index, degree)
    end
end

function nbsplines(;nknots, degree)::Int
    nknots - degree - 1
end

function (b::BSpline)(t)::Float64
    j = b.index
    d = b.degree
    ts = b.knots
    tj = ts[j]
    tj1 = ts[j+1]
    if d === 0
        return tj <= t < tj1
    end
    tjd = ts[j+d]
    tj1d = ts[j+1+d]
    return_s1 = (tj < tjd) && (tj1 == tj1d)
    return_s2 = (tj == tjd) && (tj1 < tj1d)
    if !return_s1
        s2 = (tj1d - t) / (tj1d - tj1) * BSpline(knots, j+1, d-1)(t)
    end
    if !return_s2
        s1 = (t - tj) / (tjd - tj) * BSpline(knots, j, d-1)(t)
    end
    if return_s1
        return s1
    end
    if return_s2 
        return s2
    end
    return s1 + s2
end

function MC.plottype(::BSpline)
    MC.Lines
end

function MC.convert_arguments(trait::MC.PointBased, b::BSpline)
    t_min, t_max = extrema(b.knots)
    ts = range(t_min, t_max, length=1000)
    ys = map(b, ts)
    MC.convert_arguments(trait, ts, ys)
end

d = 3
knots = [-1, 0, 1, 2, 3, 4, 5, 6]
fap = vlines(knots, color=:black)
plot!(BSpline(knots, 1, 0))
plot!(BSpline(knots, 1, 1))
plot!(BSpline(knots, 1, 2))
plot!(BSpline(knots, 1, 3))
plot!(BSpline(knots, 1, 4))
plot!(BSpline(knots, 1, 5))
# plot!(BSpline(knots, 4, 2))
fap
