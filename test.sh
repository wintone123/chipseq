awk 'BEGIN{
        p=1
    }
    {
        if(p==1||p==2){
            gsub(/^[@]/,">");print
        };
        if(p==4)p=0;p++
    }'