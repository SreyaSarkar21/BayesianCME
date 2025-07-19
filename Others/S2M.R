S2M=function(Beta,lower,upper,l)
{
    
    N=dim(Beta)[1]
    p=dim(Beta)[2]
    b.i=seq(lower,upper,length=l)
    H.b.i=NULL
    for(r in 1:l)
    {
        KK=NULL
        for (i in 1:N)
        {
            fit=kmeans(abs(Beta[i,]),2)
            cen1=min(fit$centers)
            cen2=max(fit$centers)
            temp1=Beta[i,]
            
            
            fit1=fit
            
            while(cen2-cen1>b.i[r])
            {
                fit1=fit
                temp=which.min(fit$centers)
                temp1=temp1[which(fit$cluster==temp)]
                if(length(temp1)<=2){break;}
                fit=kmeans(abs(temp1),2)
                cen2=max(fit$centers)
                cen1=min(fit$centers)
            }
            temp=which.min(fit1$centers)
            KK[i]=p-length(which(fit1$cluster==temp))
        }
        H.b.i[r]=as.numeric(names(sort(-table(KK)))[1])
    }
    abs.post.median=NULL
    for(i in 1:p)
    {
        abs.post.median[i]=median(abs(Beta[,i]))
    }
    return(list(N=N,
                p=p,
                H.b.i=H.b.i,
                b.i=b.i,
                abs.post.median=abs.post.median
    ))
}


Hbi.vs.bi=function(S2M)
{
    plot(S2M$b.i,S2M$H.b.i)
}


S2M.vs=function(S2M,H)
{
    
    p=S2M$p
    abs.post.median=S2M$abs.post.median
    var.selected=order(abs.post.median)[p:(p-H+1)]
    return(var.selected)
}