buy = 0;
buyvalue = 0;
sellvalue = 0;
sell = 1;
capital = 2000;
wallet = capital;
j = 1;
returns = 0;
% n = 7569; AAPL Full Dataset
n = 1011; % AAPL Shortened Dataset

for i=1:size(Combined,1)
    signal = Combined(i);
    if signal > 0
        if buy == 0
            buyvalue = close(i);
            buy = 1;
            sell = 0;
            wallet = wallet - buyvalue;
        end
    elseif signal < 0
        if sell == 0
            sellvalue = close(i);
            buy = 0;
            sell = 1;
            wallet = wallet + sellvalue;
            returns(j) = sellvalue - buyvalue;
            j = j+1;
        end
    end
end
%wallet
%returns = returns';
sharpe = mean(returns)/std(returns)
annual_return = ((wallet /capital)^(1/(n/250))-1) * 100
    
            
        
    