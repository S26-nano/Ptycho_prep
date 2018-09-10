function probe = v2_quick_probe(film_thickness,theta,twotheta,defocus,EkeV_in,focdetdist_in,theoryflag,rockphase,strainphase)

if(nargin<8)
    rockphase=0;strainphase=0;
end
if(nargin<6)
    theoryflag = 0; %Use experimental data, not zpdp theoretical optic
end

if(nargin<5)
    d2_bragg_in = 0; %Use same energy that probe was measured at
else
    lambda_in = 1.239842/(EkeV_in*1000); %wavelength microns
    d2_bragg_in = focdetdist_in * lambda_in/(256*55); %pixel size and range set for MPX3
end

if(nargin<4)
    defocus = 0;
end

load('/CNMshare/savedata/2017R1/Ptychography_test/Probe_library/FZP_150um_20nm/FZP_150_20_defocus_30um_9keV_mpx3.mat', 'probe_g','d2_bragg','lambda','probe_go');
defocus_probe = 0;

% If at different energy, rescale defocus value
if(d2_bragg_in>0) 
    defocus_probe = defocus_probe * lambda/(1.239842/(EkeV_in*1000));
    %display(defocus_probe)
end

% If theoretical focusing desired use zpdp calculated answer
if(theoryflag) 
    probe_g = probe_go; 
end 
    
probe = zeros(size(probe_g));

%Propagate probe from library defocus to desired defocus condition
if(defocus ~=defocus_probe)
    prbxvals(1:size(probe_g,1)) = ((1:size(probe_g,1))-1)*d2_bragg*1e-6;prbx = repmat(prbxvals,size(probe_g,2),1);
    prbyvals(1:size(probe_g,2)) = ((1:size(probe_g,2))-1)*d2_bragg*1e-6;prby = repmat(prbyvals',1,size(probe_g,1));
    [probe_g2 x_there y_there dx_there dy_there] = free_propagate_paraxial_2D(probe_g,prby,prbx,(defocus-defocus_probe)*1e-6,(2*pi/(lambda*1e-6)),(lambda*1e-6),2);
else
    probe_g2 = probe_g;
end

%Interpolate from library sampling to desired sampling - i.e. different detector distance / beam energy / pixel size
if(d2_bragg_in ~= d2_bragg)
    prbnum = size(probe_g2,2);
    prbcen = prbnum/2;
    iiq = repmat((((1:prbnum)-prbcen)*d2_bragg_in/d2_bragg+prbcen)',[1,prbnum]);
    jjq = repmat((((1:prbnum)-prbcen)*d2_bragg_in/d2_bragg+prbcen),[prbnum,1]);
    probe_g3 = interp2(probe_g2,iiq,jjq);
else
    probe_g3 = probe_g2;
end

%Project probe through a perfect thin film at desired detector angle, adding strain or rotation offsets
if(twotheta ~= 0 && film_thickness ~= 0)
    
    numpix = round((film_thickness*sind(twotheta)/sind(twotheta-theta))/d2_bragg);
    numcen = round(numpix/2);
    numprb = size(probe_g3,2);
 
    for jj = 1:numprb
        rph1 = rockphase*(jj-1)*d2_bragg/sind(twotheta-theta);
        sph1 = strainphase*(jj-1)*d2_bragg*cosd(twotheta-theta);
        kst = max([1,(jj-numcen)]);
        kend = min([numprb,(jj-numcen+numpix)]);
        for kk=kst:kend
            rph2 = rockphase*((kk-1)*d2_bragg/sind(twotheta))*cosd(twotheta-theta);
            sph2 = strainphase*((kk-1)*d2_bragg/sind(twotheta))*cosd(twotheta-theta);
            probe(:,jj)=probe(:,jj)+probe_g3(:,kk)*exp(1i*(rph1+rph2))*exp(1i*(sph1-sph2))/(kend-kst);
        end
    end
else
    probe = probe_g3;
end
%%{
 %patch for using pixirad vs mpx3
probe2 = probe(3:514,3:514);
probe = probe2;
%}        
%save('probe.mat', probe);
probe=probe(129:384,129:384);
%probe=probe(2*(1:256), 2*(1:256));
end
