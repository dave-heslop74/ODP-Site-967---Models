    
    function utils_quantile(x,q)
        xs = sort(x)
        ys = range(0,stop=1,length=length(x))
        test = LinearInterpolation(ys,xs)
    
        return test(q)
    end

    function utils_output(path,output)
        writedlm(string(path,"Fixed forcing",rand(1:100000,1)[1],".out"),  output, ',')
        return 1
    end

    function utils_quadfit(x,y)
        X = hcat(ones(length(x)),x,x.^2)
        return X\y
    end

    function utils_quadfit0(x,y)
        X = hcat(x,x.^2)
        return X\y
    end
    
    function utils_quadval(p,x)
        X = hcat(ones(length(x)),x,x.^2)
        return X*p
    end

    function utils_quadval0(p,x)
        X = hcat(x,x.^2)
        return X*p
    end

    function utils_linfit(x,y)
        X = hcat(ones(length(x)),x)
        return X\y
    end
    
    function utils_linval(p,x)
        X = hcat(ones(length(x)),x)
        return X*p
    end

    function utils_sortrows(A,dim=1)
        return A[sortperm(A[:, dim]), :]
    end

    function utils_gsmooth(X,age,sigma)
        
        output = zeros(length(age));
        for i=1:length(age)
            W = pdf.(Normal(age[i],sigma),X[:,1])
            output[i] = sum(W.*X[:,2])./sum(W)
         end
        
        return output
    end

    function utils_quadsolve(pp,y)
        
        c,b,a = pp;
        c = c.-y
        output = (-b.+sqrt.(b.^2.0.-4.0.*a.*c))/(2.0.*a);
        output = [output (-b.-sqrt.(b.^2.0.-4.0.*a.*c))/(2.0.*a)];

        return output
    end

    function utils_zscore(x)
        
        x .-= mean(x)
        x ./= std(x)

        return x
    end


    function utils_invinterp(X,Y,Y0)
        
        A = Y[1:end-1].-Y0;
        B = Y[2:end].-Y0;
        
        idx = findall(<=(0), A.*B);
        output = zeros(length(idx));

        for i=1:length(output)
            x = Y[idx[i]:idx[i]+1];
            y = X[idx[i]:idx[i]+1];
            p = sortperm(x)
            f = LinearInterpolation(x[p],y[p])
            output[i] = f(Y0);
        end
    
        return output
    
        end