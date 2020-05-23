using StatsBase;
using Logging;
using Juno;
using Random;
include("totalisticAutom.jl")
using .Totalistic
ENV["JULIA_DEBUG"] = "all"
logger = ConsoleLogger(stdout, Logging.Debug);

logspace(start,stop,length)= exp10.(range(start; stop=stop,length=length));
ncolorAxes=4;
nrules= 150;
ncells= 160;
pre_evol= 100;
ngenerations= 150;

rng= MersenneTwister(1);
rulemax= BigInt(2^ncolorAxes)^(3*2^ncolorAxes-2);
ruleset= rand(rng,big.(1:rulemax), nrules);

evolveCA(ca::TCA,gens::Int)= gens>0 ? evolveCA(nextgeneration(ca),gens-1) : ca;
produceCA(ca::TCA,gens::Int)= produceCA(nextgeneration(ca),gens,Vector{BitArray}())
produceCA(ca::TCA,gens::Int,cellList::Vector{BitArray})= gens>0 ?
    produceCA(nextgeneration(ca),gens-1, push!(cellList,ca.cells)) : cellList;

@time calist= Dict(r=>evolveCA(TCA(r,ncolorAxes,ncells),pre_evol) for r in ruleset);
@time cacells= Dict(r=>produceCA(calist[r],ngenerations) for r in ruleset);

using Luxor;
drawCA(r::BigInt)= @png begin
    sidelength=4;
    Luxor.translate(boxtopcenter(BoundingBox()) + sidelength);
    for i in 1:ngenerations
      colorDraw(cacells[r][i][:,:],sidelength,ColorSchemes.ice);
      Luxor.translate(Point(0, sidelength));
    end
  end;

## Filter out boring CAs
using Statistics, Plots;
### Autocorrelations for each rule
ac= Dict(r=>(map(g->autocor(cacells[r][g],ncells÷10:ncells÷2),1:length(cacells[r]))) for r in ruleset);
acVar= Dict(r=> mean(map(g-> log(std(ac[r][g])),1:length(ac[r])))
            for r in ruleset);
vals= d-> collect(values(sort(d)));

### Filter
filter!(p->!isnan(p.second), acVar);
ac_quant_0595= quantile(vals(acVar),[0.05,0.95]);
filter!(p-> p.second > ac_quant_0595[1] && p.second < ac_quant_0595[2], acVar);
if ncolorAxes>2
  filter!(p-> p.second > -2.38 && p.second < -2.37, acVar);
else
  filter!(p-> p.second > -2.8 && p.second < -1.9, acVar);
end
display(histogram(vals(acVar); bins=20))

using DataStructures;
acVar= OrderedDict(sort(collect(acVar), lt= (a,b)->isless(a.second,b.second)));
interestingRuleset= collect(keys(acVar));

# Display interesting CAs
using Printf;
@printf("%d interesting CA rules:\n", length(interestingRuleset));
display(interestingRuleset);
for r in interestingRuleset; drawCA(r); sleep(0.2); end
