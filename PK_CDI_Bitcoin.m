%% updates 1_9_5
% 1.EPK CDI
% 2.Rookley's Q
% 3.x-axis is Moneyness=K/S_T
%% load data
clear,clc
% addpath("m_Files_Color/colormap/")
daily_price = readtable("data/BTC_USD_yf_2017-2022.csv");

% Read the two tables
option = readtable("data/btc_option_2022.csv");

% option1 = readtable("data/final_btc_option_2022-10-01-2022-12-30.csv");
% option2 = readtable("data/final_btc_option_2022-07-01-2022-09-29.csv");
% option = vertcat(option1, option2);


%% define tau 
unique_tau_values = unique(option.tau);
tau_values = unique_tau_values(unique_tau_values >= 7 & unique_tau_values <= 60)';
tau_values = [7,8,13,30]; %

%% for loop for all tau values
% Loop over each tau value
for tau = tau_values
    % dates from Q density files
    files = dir(fullfile('Q_density_files_new', sprintf('*_tau%d.csv', tau)));
    dates = strings(1, length(files));
    for i = 1:length(files)
        dateStr = regexp(files(i).name, '\d{4}-\d{2}-\d{2}', 'match', 'once');
        if ~isempty(dateStr)
            dates(i) = dateStr;
        end
    end
    dates_Q = dates';
    
    % Dates from option data
    dates_option = string(datetime(datestr(unique(option.date((option.tau >= tau) & (option.tau <= tau))), 'YYYY-mm-DD'), "Format", "yyyy-MM-dd"));
    dates = intersect(dates_option, dates_Q);
    
    dates_CDI = dates(2:end);
    dates_Q = dates(1:(end-1));
    
    %% Estimate EPK
    figure;
    for i0=1:numel(dates_CDI)
        [~,~,~]=mkdir(strcat("EPK_figures/"));
        
        % prepare Q-density
        % dates_Q_string = string(datetime(dates_Q,"Format","yyyyMMdd"));
        % dates_cell = dates_Q_string(1:i0);
        dates_cell = dates_Q(1:i0);
    
        realizedKhRet=cell(length(dates_cell),1);
        realizedQdenRet=cell(length(dates_cell),1);
    
        for i = 1:length(dates_cell)
    
            a = sprintf("Q_density_files_new/raw_Q_density_%s_tau%d.csv", dates_cell(i), tau);
            data_q = readtable(a);
    
            sp1=daily_price;
            sp1(datenum(sp1.Date)<datenum(dates_cell(i),"yyyy-MM-dd") | datenum(sp1.Date)>datenum(dates_cell(i),"yyyy-MM-dd")+tau,:)=[];
            rt=linspace(-1,log(sp1.Adj_Close(end)/sp1.Adj_Close(1)),500);
    
            Qdentisty = spline(data_q.m,data_q.y,rt);
    
            realizedKhRet{i,1}=rt;
            realizedQdenRet{i,1}=Qdentisty;
        end
    
        % delete empty cell if any
        i = 1;
        while i <= length(realizedKhRet)
            if isempty(realizedQdenRet{i,1})
                realizedQdenRet(i,:)=[];
                realizedKhRet(i,:)=[];
            else
                i=i+1;
            end
        end
        
        % CDI part 4, 4, 0 .0001 end of day price
        [sampleestimate_4_4_00001, returns_4_4_00001] = CDI_estimator(realizedKhRet,realizedQdenRet, @OptSDF,4,4, .0001) ;
    
        % plot the moneyness
        Colors = rainbow(7);
        ret=returns_4_4_00001 ;
        moneyness = exp(ret);
        plot(moneyness, sampleestimate_4_4_00001,"Color", 'k', 'LineWidth', 2)
    
        hold on
        [sampleestimate_5_5_00001, returns_5_5_00001] = CDI_estimator(realizedKhRet,realizedQdenRet, @OptSDF,5,5, .0001);
        ret=returns_5_5_00001;
        moneyness = exp(ret);
        plot(moneyness, sampleestimate_5_5_00001,"Color",Colors(2,:), 'LineWidth', 2);
        [sampleestimate_6_6_00001, returns_6_6_00001] = CDI_estimator(realizedKhRet,realizedQdenRet, @OptSDF,6,6, .0001);
        ret=returns_6_6_00001;
        moneyness = exp(ret);
        plot(moneyness, sampleestimate_6_6_00001,"Color",Colors(3,:), 'LineWidth', 2);
        [sampleestimate_7_7_00001, returns_7_7_00001] = CDI_estimator(realizedKhRet,realizedQdenRet, @OptSDF,7,7, .0001);
        ret=returns_7_7_00001;
        moneyness = exp(ret);
        plot(moneyness, sampleestimate_7_7_00001,"Color",Colors(4,:), 'LineWidth', 2);
    
        xlabel('Moneyness'),ylabel('EPK')
        xlim([0.75,1.02]);
        % ylim([0.5 3]);     
        hold off
        % legend('4, 4','5, 5','6, 6','7, 7');

        % title(sprintf('%s, tau = %d, # moments, # knots', dates_CDI(i0), tau))
        saveas(gcf, sprintf("EPK_figures/BRC_EPK_%dday_%s.png", tau, dates_CDI(i0)))

        % outtable = table(ret',sampleestimate_4_4_00001,sampleestimate_5_5_00001,sampleestimate_6_6_00001, ...
        %     sampleestimate_7_7_00001,sampleestimate_5_8_00001,sampleestimate_5_9_00001,sampleestimate_5_10_00001, ...
        %     'VariableNames',{'Return','J_4_m_4','J_5_m_5','J_6_m_6','J_7_m_7','J_5_m_8','J_5_m_9','J_5_m_10'});
        % writetable(outtable,strcat("EPK_figures/CDI_Rookley_1_9_0_return/EPK_9day_",dates_CDI(i0),".csv"))
    end
end

