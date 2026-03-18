function [elecmatrix,eleclabels,anatomy]=getOPSCEAelecs(pt,clin1TDT2MNI3)

	%data_root = getenv("KLEEN_DATA");
    data_root='/data/';
	mainpath=fullfile(data_root, 'opscea');      %path for OPSCEA folders
	% if exist(mainpath,'dir')==0
	% 	mainpath='/Volumes/KLEEN_DRIVE/OPSCEA'; 
	% end
	% if ~(exist(mainpath,'dir')==7) %jk put parenthesis to make it a falsity of entire statement
	% 	msgbox('KLEEN_DRIVE not found, you may need to reconnect?');
	% 	return;
	% end
	
	ptdir=fullfile(mainpath, 'OPSCEADATA', pt);
	filepath='/Imaging/elecs/';
	if clin1TDT2MNI3==1
		typ='clinical';
	elseif clin1TDT2MNI3==2
		typ='TDT';
	elseif clin1TDT2MNI3>=3
		typ='TDT';
	end
	clinfilename_alternate='clinical_TDT'; 
	ext='_elecs_all.mat'; 
	elecfile = 'Electrodefile.mat';
	
	
	if clin1TDT2MNI3<3
		fn=fullfile(ptdir, filepath, [typ ext]);
	elseif clin1TDT2MNI3==3
		fn=fullfile(ptdir, filepath, [typ '_elecs_all_warped.mat']);
	elseif clin1TDT2MNI3==4
		fn=fullfile(ptdir, filepath, [typ '_elecs_all_yale.mat']);
	end
	fn_alt=fullfile(ptdir, filepath, clinfilename_alternate, ext);
	fn_elec=fullfile(ptdir, filepath, elecfile);
	
	if exist(fn)==2     
		load(fn); 
	elseif exist(fn_alt)==2
		load(fn_alt); 
	elseif exist(fn_elec)==2
		load(fn_elec);
	else elecmatrix=[];
		eleclabels=[];
		anatomy=[];
		disp('No result, check path and filename');
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