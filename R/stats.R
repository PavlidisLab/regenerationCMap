#' @export
a = function(V,n){
    t = length(V)
    1:t %>% sapply(function(j){
        g=j/t - V[j]/n
    }) %>% max
}

#' @export
b = function(V,n){
    t = length(V)
    1:t %>% sapply(function(j){
        V[j]/n - (j-1)/t
    }) %>% max
}


#' @export
ks = function(a,b){
    if(a>b){
        return(a)
    } else if (b>a){
        return(-b)
    } else if(a == b){
        return(0)
    }
}

ksCalc = function(V,n){
    a = a(V,n)
    b = b(V,n)
    return(ks(a,b))
}

#' @export
score = function(kUp,kDown){
    if(sign(kUp) == sign(kDown)){
        return(0)
    } else {
        kUp - kDown
    }
}

#' @export
scoreCalc = function(Vup,Vdown,n){
    aUp = a(Vup,n)
    bUp = b(Vup,n)
    kUp = ks(aUp,bUp)
    
    aDown = a(Vdown,n)
    bDown = b(Vdown,n)
    kDown = ks(aDown,bDown)
    score = score(kUp,kDown)
    return(c(kUp = kUp, kDown = kDown,score = score))
}
