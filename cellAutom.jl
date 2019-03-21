# Cellular automata
# Source: https://github.com/cormullion/cormullion.github.io/blob/master/source/automata.jl

module SimulatorCA

mutable struct CA
  rule::Int64
  cells::BitArray{1}
  colorstops::Array{Float64, 1}
  ruleset::BitArray{1}
  generation::Int64
  function CA(rule, ncells = 100)
    cells                    = falses(ncells)
    colorstops               = zeros(Float64, ncells)
    ruleset                  = binary_to_array(rule)
    cells[length(cells) ÷ 2] = true
    generation               = 1
    new(rule, cells, colorstops, ruleset, generation)
  end
end
export CA;

function binary_to_array(n)
  a = BitArray{1}()
  for c in 7:-1:0
    k = n >> c
    push!(a, (k & 1 == 1 ? true : false))
  end
  return a
end
function rules(ca::CA, a, b, c)
  lng = length(ca.ruleset)
  return ca.ruleset[mod1(lng - (4a + 2b + c), lng)]
end
function nextgeneration(ca::CA)
  l = length(ca.cells)
  nextgen = falses(l)
  for i in 1:l
    left   = ca.cells[mod1(i - 1, l)]
    me     = ca.cells[mod1(i, l)]
    right  = ca.cells[mod1(i + 1, l)]
    nextgen[i] = rules(ca, left, me, right)
  end
  ca.cells = nextgen
  ca.generation += 1
  return ca
end
export nextgeneration;

Base.show(io::IO, ::MIME"text/plain", ca::CA) =
print(io, join([c ? "■" : " " for c in ca.cells]))


using Luxor, Colors

function draw(cells::BitArray, cellwidth=10)
  lng = length(cells)
  for i in 1:lng
    if cells[i] == true
      pt = Point(-(lng ÷ 2) * cellwidth + i * cellwidth, 0)
      box(pt, cellwidth, cellwidth, :fill)
    end
  end
end
export draw;

function drawManyRules(ruleset::AbstractArray)
  for rule_i in ruleset
    display(rule_i)
    @png begin
      ca = CA(rule_i, 200)
      sidelength = 4
      ## start at the top
      translate(boxtopcenter(BoundingBox()) + sidelength)
      for _ in 1:200
        nextgeneration(ca);
      end
      for i in 1:200
        draw(ca.cells, sidelength)
        nextgeneration(ca)
        translate(Point(0, sidelength))
      end
    end #800 850 "images/automata/simple-ca.png"
    sleep(0.15);
  end
end
# rule45

end #module
