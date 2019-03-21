using StatsBase;
using Main.SimulatorCA;

ruleset=1:255;
evolveCA(ca::CA,gens::Int)= gens>0 ? evolveCA(nextgeneration(ca),gens-1) : ca;
produceCA(ca::CA,gens::Int)= produceCA(nextgeneration(ca),gens,Vector{BitArray}())
produceCA(ca::CA,gens::Int,cellList::Vector{BitArray})= gens>0 ?
    produceCA(nextgeneration(ca),gens-1, push!(cellList,ca.cells)) : cellList;
calist= Dict(r=>evolveCA(CA(r,200),200) for r in ruleset);
cacells= Dict(r=>produceCA(calist[r],200) for r in ruleset)

using Luxor;
drawCA(r::Int)= @png begin
    sidelength=4;
    Luxor.translate(boxtopcenter(BoundingBox()) + sidelength);
    for i in 1:200
      draw(cacells[r][i],sidelength);
      Luxor.translate(Point(0, sidelength));
    end
  end;

## Filter out boring CAs
using Statistics, Plots;
ac= Dict(r=>(map(g->autocor(cacells[r][g],30:100),1:length(cacells[r]))) for r in ruleset);
acVar= Dict(r=> mean(map(g-> log(std(ac[r][g])),1:length(ac[r])))
            for r in ruleset);
filter!(p->!isnan(p.second), acVar);
vals= d-> collect(values(sort(d)));
ac_quant_0595= quantile(vals(acVar),[0.05,0.95]);
filter!(p-> p.second > ac_quant_0595[1] && p.second < ac_quant_0595[2], acVar);
filter!(p-> p.second > -3 && p.second < -1.8, acVar);
display(histogram(vals(acVar); bins=20))

using DataStructures;
acVar= OrderedDict(sort(collect(acVar), lt= (a,b)->isless(a.second,b.second)));
interestingRuleset= collect(keys(acVar));

# Display interesting CAs
using Printf;
@printf("%d interesting CA rules:\n", length(interestingRuleset));
display(interestingRuleset);
for r in interestingRuleset; drawCA(r); sleep(0.2); end
