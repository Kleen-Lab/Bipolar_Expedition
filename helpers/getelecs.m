function [elecmatrix,eleclabels,anatomy]=getelecs(pt,clin1TDT2MNI3)

do_opscea_check = false;

%data_root = '/Volumes/SPIKE/data/';
%if strcmp(data_root, '')
%	disp("It seems that the $KLEEN_DATA environment variable is not set, you might want to set it and mount the data.");
%end

data_root = '/data/';

mainpath=fullfile(data_root, 'imaging');
% if exist(mainpath,'dir')==0
% 	mainpath='/Volumes/KLEEN_DRIVE/AN/DATA'; 
% end
% if ~(exist(mainpath,'dir')==7) %jk put parenthesis to make it a falsity of entire statement
% 	msgbox('KLEEN_DRIVE not found, you may need to reconnect?');
% 	return;
% end


%disp(mainpath);

mainpath='/data/imaging'

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

% one last check, in the Chang Lab server
if isempty(elecmatrix);
    %disp('checking Chang Lab server, if connected')
    CL_filepath=fullfile('~/DataMount/data_store2/imaging/subjects/', pt, 'elecs', [typ ext]);
    CL_filepath_alt=fullfile('~/DataMount/data_store2/imaging/subjects/', pt, 'elecs', [clinfilename_alternate ext]);
    if exist(CL_filepath)==2
		load(CL_filepath);
        %disp(CL_filepath)
    elseif exist(CL_filepath_alt)==2
		load(CL_filepath_alt);
        %disp(CL_filepath_alt)
	end
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