
function Tlist = count2list(n)

Tlist=[];
for i=1:max(n)
    Tlist = [Tlist; find(n>=i)];
end
Tlist =sort(Tlist);