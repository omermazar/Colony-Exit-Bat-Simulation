fn1 = fieldnames(Params);
ix=1;
mm=0;
for n = 1:numel(fn1)
    cn = fn1{n};
    fn2 = fieldnames(Params.(cn));
    for m = 1:numel(fn2)
        cnn = fn2{m};
        if isnumeric(Params.(cn).(cnn))
            if numel(Params.(cn).(cnn)) == 1
                t(ix) = Params.(cn).(cnn) == AllParams.(cn).(cnn);
            else % numel
                t(ix) = false;
            end % numel
        else % if
            t(ix) = strcmp(Params.(cn).(cnn) , AllParams.(cn).(cnn) );
        end % if
        if ~t(ix)
            mm = mm+1;
            err.sub(mm) = string(cnn);
            err.main(mm) = string(cn);
        end
        ix= ix+1;
    end % for m
end %for n 