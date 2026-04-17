function [sigcolor,nonsigcolor] = generate_distinguishablecolors(color1,color2,numcolors)
    if ~isa(color1,'double')
    color1 = hex2rgb(color1);
    color2 = hex2rgb(color2);
    end
    orderededges = [color1;color2];
    orderedcolors = zeros(numcolors+1,size(color1,2));
    orderedcolors(1:2,:) = orderededges;
    currnumcolors = 2;
    while currnumcolors<(numcolors+1)
        newedges = zeros(2*size(orderededges,1)-1,size(orderededges,2));
        newedges(1:2:end,:) = orderededges;
        edgecnt = 0;
        for jj =1:size(orderededges,1)-1
            edgecnt = edgecnt+2;
            curredges = orderededges(jj:jj+1,:);
            currnewedge = mean(curredges);
            newedges(edgecnt,:) = currnewedge;
            try
                currnumcolors = currnumcolors+1;
                orderedcolors(currnumcolors,:) = currnewedge;
            catch
                break
            end
        end
        orderededges = newedges;
    end
    sigcolor = cellstr(rgb2hex(orderedcolors(1:numcolors,:)));
    nonsigcolor = cellstr(rgb2hex(orderedcolors(numcolors+1,:)));
end