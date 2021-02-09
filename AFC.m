function [assimilate,melts,Composition,Mass,Conditions]=AFC(bulk,As_Names,As_Comp,As_Prop,As_Cond,X,len,dh,R)
% AFC function to calculate the chemical consequences of assimilation of
% magmatic mush. This code uses rhyoliteMELTS v1.2.0 (Gualda and Ghiorso,
% 2015). Please ensure that you have added the relevant alphaMELTS for
% MATLAB folders to your path. This code uses the alphaMELTS for MATLAB
% functions as of Dec 2020.
%   Required Inputs: 
%           'bulk' (19-by-1 vector with initial liquid composition)
%           'As_Names' (string of length n with the names you assign to the
%           assimilated minerals; e.g. "olivine_as")
%           'As_Comp' (19-by-n matrix with the composition of each mineral 
%           in the assimilated component; where n is the number of
%           minerals)
%           present in the assimilated component)
%           'As_Prop' (1-by-n vector with the proportions of each mineral
%           in the assimilated component. This does not need to equal 1 or
%           100)
%           'As_Cond' (2-by-1 vector containing pressure (1) and
%           temperature (2))
%           'X' (g of solid to assimilate at each step (per 1g of melt))
%           'len' (number of loops)
%   Optional Inputs:
%           'dh' (change in enthalpy of the system (per 100g) at each step)
%
%   Outputs:
%           'assimilate' (MELTSdynamic for the assimilated components)
%           'melts' (MELTSdynamic for the assimilation path)
%           'Composition' (Stacked table containing the composition of each
%           phase that is present in the system; for example,
%           Composition.liquid1.MgO would return the MgO content of the
%           liquid at each step)
%           'Mass' (Table containing the mass (per 100g) of each phase at
%           each step)
%           'Conditions' (Table containing the pressure and temperature at
%           each model step)
%       Note: If only assimilate and melts are specified as outputs the
%       code will terminate after the main for loop.

    if nargin < 8 
        dh=0;
    end
    if nargin < 9
        R=0;
    end

    %chose MELTS dynamic
    assimilate = MELTSdynamic(1);
    
    % mass of assimilant at each step (g added to 1g of melt)
    X_a=X;
    
    % conditions of assimilant
    assimilate.engine.pressure=As_Cond(1);
    assimilate.engine.temperature=As_Cond(2);
    
    % calculate phase properties
    assimilate.engine.calcPhaseProperties(As_Names,As_Comp);
    
    H_As=zeros(1,length(As_Names));
    X_As=zeros(19,length(As_Names));
    for i=1:length(As_Names)
        H_As(i)=assimilate.engine.getProperty('h',As_Names(i));
        X_As(:,i)=assimilate.engine.getProperty('dispComposition',As_Names(i));        
    end
    
    % define composition and enthalpy of assimilant
    H_As_f=sum(As_Prop.*H_As)./sum(As_Prop);
    X_As_f=sum(As_Prop.*X_As,2)./sum(As_Prop);
    
    % ensure enthalpy and composition are relative to 100g
    H_As_f=100.*H_As_f./sum(X_As_f);
    X_As_f=100.*X_As_f./sum(X_As_f);
    
    % normalise bulk composition to 100g
    liq=bulk/sum(bulk)*100;
    
    % define new melts dynamic
    melts=MELTSdynamic(1);
    melts.engine.pressure=As_Cond(1);
    melts.engine.temperature=As_Cond(2);

    % set bulk composition as initial liquid
    melts.engine.setBulkComposition(liq);

    % find liquidus
    melts.engine.findLiquidus;
    
    % calc equilibrium
    melts.engine.calcEquilibriumState(1,0);

    % calc enthaly of system
    Melt_H=melts.engine.h('bulk');
    Melt_H_In=Melt_H;
    
    % loop to simulate the influence of progressive reactive porous flow
    for i=1:len
        % add new node
        melts=melts.addNodeAfter;
        
        % calculate enthalpy at this step based on the assimilant addition
        H=dh+(Melt_H+H_As_f*X_a)/(1+X_a); % enthalpy per 100 g
        
        % set new bulk composition (ensure it's 100g)
        newBulk=liq+X_As_f*X_a;
        newBulk=newBulk/sum(newBulk)*100;
        
        % set new bulk composition of the system
        melts.engine.setBulkComposition(newBulk);
        
        % calc equilibrium state
        melts.engine.calcEquilibriumState(1,0);
        
        disp(melts.engine.status.message);
        tf = isequal(melts.engine.status.message,'Quadratic iterations exceeded.');
        if tf ==1
            break
        end
        % calc enthalpy (at this moment the temperature is equivalent to the previous step)
        H1=melts.engine.h("bulk");
        
        % iterate temperature to find solution at constant enthalpy
        if H1<H
            while H1<H
                melts.engine.temperature=melts.engine.temperature+0.2;
                melts.engine.calcEquilibriumState(1,0);
                H1=melts.engine.h("bulk");
            end
        end

        if H1>H
            while H1>H
                melts.engine.temperature=melts.engine.temperature-0.2;
                melts.engine.calcEquilibriumState(1,0);
                H1=melts.engine.h("bulk");
            end
        end
        
        % extract new melt composition
        liq=melts.engine.getProperty('dispComposition','liquid1');
        liq=liq/(sum(liq))*100;
        
        % calculate composiiton of melt after re-injection
        liq=(1-R).*liq+R.*(bulk/sum(bulk)*100);
        
        if liq(19)~=0
            break
        end

        % determine entropy of the liquid phase and set to enthalpy per 100g
        Melt_H=melts.engine.h("liquid1")/melts.engine.getProperty('mass','liquid1')*100;
        Melt_H=Melt_H.*(1-R)+Melt_H_In.*R;
    end
    
    if nargout>2
        len=length(melts.getListProperty('pressure'))-1;
        % extract conditions of model
        N={'P','T'};
        Conditions=table(zeros(len+1,1),zeros(len+1,1));
        Conditions.Properties.VariableNames=N;
        Conditions.P=melts.getListProperty('pressure')';
        Conditions.T=melts.getListProperty('temperature')';

        % Define names of all possible phases in the MELTS models (this can be
        % cut down to save time
        Names={'liquid1','fluid1','olivine1','olivine2','clinopyroxene1',...
            'clinopyroxene2','plagioclase1','plagioclase2',...
            'spinel1','spinel2','orthopyroxene1','orthopyroxene2',...
            'kfeldspar1','kfeldspar2','apatite1','rhmoxide1',...
            'quartz1','biotite1','whitlockite1'};
        Mass=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
        Mass.Properties.VariableNames=Names;

        SolidPhases=zeros(len+1,19,width(Mass));
        Elements={'SiO2','TiO2','Al2O3','Fe2O3','Cr2O3','FeO','MnO','MgO',...
            'NiO','CoO','CaO','Na2O','K2O','P2O5','H2O','CO2','SO2','Cl','F'};
        Phase=table(zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),...
            zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),...
            zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),...
            zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),...
            zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),zeros(len+1,1));
        Phase.Properties.VariableNames=Elements;
        Composition=table(Phase,Phase,Phase,Phase,...
            Phase,Phase,Phase,Phase,Phase,...
            Phase,Phase,Phase,Phase,Phase,...
            Phase,Phase,Phase,Phase,Phase);
        Composition.Properties.VariableNames=Names;

        for j=1:width(Mass)
            Mass{1:len+1,j}=melts.getListProperty('mass',Names(j))';
            SolidPhases(:,:,j)=melts.getListProperty('dispComposition',Names(j))';
            Phase{:,:}=SolidPhases(:,:,j);
            Composition(:,j)=table(Phase);
            if Names(j)=="rhmoxide1"
                Mass{1:len+1,j}=melts.getListProperty('mass','rhm-oxide1')';
                SolidPhases(:,:,j)=melts.getListProperty('dispComposition','rhm-oxide1')';
                Phase{:,:}=SolidPhases(:,:,j);
                Composition(:,j)=table(Phase);
            end 
            if Names(j)=="kfeldspar1"
                Mass{1:len+1,j}=melts.getListProperty('mass','k-feldspar1')';
                SolidPhases(:,:,j)=melts.getListProperty('dispComposition','k-feldspar1')';
                Phase{:,:}=SolidPhases(:,:,j);
                Composition(:,j)=table(Phase);
            end  
            if Names(j)=="kfeldspar2"
                Mass{1:len+1,j}=melts.getListProperty('mass','k-feldspar2')';
                SolidPhases(:,:,j)=melts.getListProperty('dispComposition','k-feldspar1')';
                Phase{:,:}=SolidPhases(:,:,j);
                Composition(:,j)=table(Phase);
            end  
        end
    end
end