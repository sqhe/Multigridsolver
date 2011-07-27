function ret = assert( b )
    if( ~b )
        if( nargin<2 || isempty(s) )
            s = ' '; 
        end
        error( s ); 
    end
end