module Totalistic

struct TCA
  rule::BigInt
  ruleset::Array{BitArray{1},1}  #ncolor-base numerical representation of the rule
  cells::BitArray{2}
  generation::Int
  ncolorAxes::Int8
  ncells::Int
  colorWeights::Array{Int,2}

  function TCA(rule::BigInt,ncolorAxes::Int,ncells::Int=100)
    cells= falses(ncells,ncolorAxes);
    cells[ncells ÷ 2,:]= trues(ncolorAxes);
    generation= 1;
    ruleset= kary_rule_transform(rule,ncolorAxes);
    colorWeights= weighColors(ncolorAxes);
    new(rule, ruleset, cells, generation, ncolorAxes, ncells, colorWeights);
  end
  TCA(cp::TCA,nextgen)= new(cp.rule,cp.ruleset,nextgen,cp.generation,
                            cp.ncolorAxes,cp.ncells,cp.colorWeights);
end
TCA(cp::TCA)= TCA(cp,cp.cells);
export TCA;

# Transform rule into an ncolor-arity number dict
function kary_rule_transform(rule::BigInt,ncolorAxes::Int)
  ncolor= 2^ncolorAxes;
  # LSB..MSB
  kary_transform(n,k)::Vector{UInt8}= begin
    kary= Vector{UInt8}();
    push!(kary, n % k);
    n= n ÷ k;
    while n>0
      push!(kary, n % k);
      n= n ÷ k;
    end
    return kary;
  end
  pad_dims(k::BitArray{1})= length(k) < ncolorAxes ?
                              [k;falses(ncolorAxes-length(k))] : k;
  colorAxesRep(k::UInt8)::BitArray{1}= pad_dims(kary_transform(k,2) .!= 0);
  colorAxesRep(kary::Vector{UInt8})::Array{BitArray{1},1}=
      [colorAxesRep(k) for k in kary];
  kary= kary_transform(rule,ncolor);
  return colorAxesRep(kary);
end

weighColors(ncolorAxes,colorWeights=[1])= ncolorAxes>1 ?
    weighColors(ncolorAxes-1,[colorWeights colorWeights[end]*2]) : colorWeights;

function rules(tca::TCA, a,b,c)::BitArray{1}
  avg_col= Int((tca.colorWeights*(a+b+c))[1]);
  nextColor= avg_col+1 > length(tca.ruleset) ?
              falses(length(a)) : tca.ruleset[avg_col+1];
  return nextColor;
end

function nextgeneration(tca::TCA)::TCA
  nextgen= falses(size(tca.cells));
  for i in 1:tca.ncells
    left   = tca.cells[mod1(i - 1, tca.ncells),:];
    me     = tca.cells[mod1(i, tca.ncells),:];
    right  = tca.cells[mod1(i + 1, tca.ncells),:];
    nextgen[i,:] = rules(tca, left, me, right);
  end
  return TCA(tca,nextgen)
end
export nextgeneration;

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

using ColorSchemes;
function colorDraw(cells::BitArray, cellwidth=10, colormap=ColorSchemes.phase)
  lng= size(cells,1); ncolorAxes= size(cells,2);
  colorWeights= weighColors(ncolorAxes);
  for i in 1:lng
    if length(colorWeights) > 1
      colorIdx= (colorWeights*cells[i,:]) ./
                  (colorWeights*ones(size(cells[i,:])));
    else
      colorIdx= cells[i] ? 1 : 0;
    end
    sethue(get(colormap, colorIdx[1]))
    pt = Point(-(lng ÷ 2) * cellwidth + (i * cellwidth), 0)
    box(pt, cellwidth, cellwidth, :fill)
  end
end
export colorDraw;

end #module
