#' @export
a =function(V,n){
    V %<>% as.matrix
    t = nrow(V)
    if(nrow(V)==0){
        return(0)
    }
    ((1:t)/t - V/n) %>% apply(2,max)
}


#' @export
b = function(V,n){
    V %<>% as.matrix
    t = nrow(V)
    if(nrow(V)==0){
        return(0)
    }
    (V/n - (1:t - 1)/t) %>% as.matrix%>% apply(2,max)
}


#' @export
ks = function(a,b){
    (a>b)*a + (b>a)*(-b)
}

#' @export
ksCalc = function(V,n){
    a = a(V,n)
    b = b(V,n)
    return(ks(a,b))
}


#' @export
score = function(kUp,kDown){
    (sign(kUp)!=sign(kDown))*(kUp-kDown)
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
    return(data.frame(kUp = kUp, kDown = kDown,score = score,instance= colnames(Vup)))
}
