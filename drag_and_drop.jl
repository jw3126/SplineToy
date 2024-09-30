using GLMakie
import MakieCore as MC
using ArgCheck
using StaticArrays

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

(;fig, ax, positions) = make_drag_and_drop()
spline = map(Hermite, positions)
lines!(ax, spline, color=:green)
fig
