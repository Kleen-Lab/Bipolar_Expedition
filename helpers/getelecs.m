function [elecmatrix,eleclabels,anatomy]=getelecs(pt,clin1TDT2MNI3)

do_opscea_check = false;

data_root = getenv("BIPOLAR_DATA");

mainpath=fullfile(data_root, 'imaging');

ptdir=fullfile(mainpath, pt);
filepath=fullfile('elecs');
if clin1TDT2MNI3==1
	typ='clinical'; 
elseif clin1TDT2MNI3==2
	typ='TDT';
elseif clin1TDT2MNI3>=3
	typ='TDT';
end

if pt == 'EC222'
    clinfilename_alternate='TDT'; 
else
    clinfilename_alternate='clinical_TDT'; 
end

if clin1TDT2MNI3 == 2
    ext = '_elecs_all.mat';
elseif clin1TDT2MNI3 == 3
    ext='_elecs_all_warped.mat'; %change for MNI brain
else
end

%disp(clin1TDT2MNI3);


if clin1TDT2MNI3<3
	fn=fullfile(ptdir, filepath, [typ ext]);
    %disp(fn)
else
	fn=fullfile(ptdir, filepath, [typ '_elecs_all_warped.mat']);
    if ~exist(fn, 'file')
        anpath = fullfile(data_root, 'an', 'preprocessed');
        ptdir = fullfile(anpath, pt);
        fn = fullfile(ptdir, filepath, [typ '_elecs_all_warped.mat']);
    end
end

if clin1TDT2MNI3==4
    fn(end-9:end)=[]; fn=[fn 'yale.mat'];
end

fn_alt=fullfile(ptdir, filepath, [clinfilename_alternate, ext]);

%disp(fn);
if exist(fn, 'file')
    load(fn);
elseif exist(fn_alt)==2
	load(fn_alt);
else
    if do_opscea_check
        [elecmatrix,eleclabels,anatomy]=getOPSCEAelecs(pt,clin1TDT2MNI3); %check OPSCEA pts
        if ~exist('elecmatrix','var')
            elecmatrix=[]; eleclabels=[]; anatomy=[];
		    disp('No result, check path and filename')
        end
    end
end

if ~exist('elecmatrix','var')
    elecmatrix=[];
end


if ~exist('eleclabels','var')
	eleclabels=[];
end
if ~exist('anatomy','var')
	anatomy=eleclabels;
	anatomy=[anatomy cell(size(anatomy,1),1)];
end
if clin1TDT2MNI3==4
    anatomy=labels; clear labels
end