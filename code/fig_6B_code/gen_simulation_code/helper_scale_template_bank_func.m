function [ scaled_template_bank ] = helper_scale_template_bank_func( template_bank, scaling_ratios )
%Funcation that scales the template voltage to certain ratios
%   Input: template_bank: 1 by n cell array with templates of hand picked
%   templates
%   Input: scaling ratio: 1 by n vector with ratios that scales the
%   template bank
%   Output: scaled_template_bank: 1 by n cell array with template_bank
%   scaled to its scaled ratio

scaled_template_bank = template_bank;

for i = 1:length(scaling_ratios)
    
    scaled_template_bank{i} = template_bank{i} * scaling_ratios(i);
        
end

end

