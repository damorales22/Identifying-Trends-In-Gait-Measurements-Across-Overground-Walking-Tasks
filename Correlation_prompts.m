function [correlation] = Correlation_prompts(var, prompt0, prompt1, prompt2)
        
    y_values = [var.(prompt0)(:)', var.(prompt1)(:)', var.(prompt2)(:)'];
    x_values = [0*ones(1,length(var.(prompt0)(:)')), 1*ones(1,length(var.(prompt1)(:)')), 2*ones(1,length(var.(prompt2)(:)'))]; % It doesn't matter what you use for the length, it's all the same

    correlation = corr(x_values(:),y_values(:));
end