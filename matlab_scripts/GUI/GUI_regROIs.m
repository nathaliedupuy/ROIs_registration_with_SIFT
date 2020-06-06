% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%           INPUTS / OUTPUT
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% GUI_regROIs(REG_DATA, REG_RES, output_var_name)
%
% REG_DATA      : struct with the TIFFs and ROIs masks
% REG_RES       : struct with the results (i.e. matched_ROIs, scores, etc)
% output_var_name   : string - name of variable to output modifications 
%
% Returns an empty array (no output), but saves the modified structure
% 'REG_RES' in a variable named 'output_var_name' using assignin. If the
% variabke does not exist, a new variable will be created in the base
% workspace; otherwise overwrites it.
%
% NOTE: to get an output from GUI, could have used this
% https://blogs.mathworks.com/videos/2010/02/12/advanced-getting-an-output-from-a-guide-gui/
% (also see MATLAB GUI example 'Modal Question Dialog')
% >> BUT cannot use UIWAIT otherwise interaction with UI tools like zoom do
% not work anymore / will use assignin instead, although not ideal solution
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%           LIST FUNCTIONS
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% GUI - MAIN
% -----------
% GUI_regROIs               : Main function - DO NOT EDIT (matlab)
% GUI_regROIs_OpeningFcn    : Init the GUI. Executes just before
%                             GUI_regROIs is made visible.
% GUI_regROIs_OutputFcn     : Return empty
% uifigure_regrois_CloseRequestFcn   
%                           : Saves the modifcations of REG_RES into
%                             specified variable in base workspace before
%                             close the figure
% uifigure_regrois_SizeChangedFcn
%
% INITIALISATION FUNCTIONS
% -------------------------
%
% SHOW/HIDE ROIS \
%   - CHOOSE TYPE TO DISPLAY: PAIRED, NON-MATCHED, OUT_OF_FOV, DELETED
%   - SHOW BOTH SESSIONS ON TIFF
% ---------------------------------------------------------------------
% toogle_rois_visibility






%__________________________________________________________________________
%
%   GUI - MAIN
%__________________________________________________________________________

function varargout = GUI_regROIs(varargin)
% GUI_REGROIS MATLAB code for GUI_regROIs.fig
%      GUI_REGROIS, by itself, creates a new GUI_REGROIS or raises the existing
%      singleton*.
%
%      H = GUI_REGROIS returns the handle to a new GUI_REGROIS or the handle to
%      the existing singleton*.
%
%      GUI_REGROIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_REGROIS.M with the given input arguments.
%
%      GUI_REGROIS('Property','Value',...) creates a new GUI_REGROIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_regROIs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_regROIs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_regROIs

% Last Modified by GUIDE v2.5 05-Apr-2020 11:33:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_regROIs_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_regROIs_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before GUI_regROIs is made visible.
function GUI_regROIs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_regROIs (see VARARGIN)

% Choose default command line output for GUI_regROIs
handles.output = hObject;

% Pass inputs into arguments
handles.REG_DATA=varargin{1};
handles.REG_RES=varargin{2};
handles.output_var=varargin{3};

% Default units should be pixels - but just in case
set(handles.uifigure_regrois,'Units','pixels')
set(handles.uiaxes1,'Units','pixels')
set(handles.uiaxes2,'Units','pixels')
set(handles.uipanel_inspect,'Units','pixels')
set(handles.uipanel_view,'Units','pixels')
set(handles.uipanel_showrois,'Units','pixels')
set(handles.uipanel_editresults,'Units','pixels')
set([handles.popupmenu_s1, handles.popupmenu_s2, handles.popupmenu_inspect, ...
    handles.listbox_rois_s1, handles.listbox_rois_s2],'Units','pixels')
set([handles.pushbutton_resetview, handles.togglebutton_hidetiff1, handles.togglebutton_hidetiff2, ...
    handles.uibuttongroup_alignmask, handles.text_linkaxes, handles.radiobutton_s1, handles.radiobutton_s2, handles.checkbox_linkaxes],'Units','pixels')
set([handles.text_distance, handles.text_distancevalue, handles.pushbutton_match, ...
    handles.pushbutton_deletepair, handles.pushbutton_deleteroi, handles.pushbutton_addroi, handles.pushbutton_outoffov],'Units','pixels')
set([handles.checkbox_showboth, handles.text_showboth, handles.checkbox_showdeleted, ...
    handles.checkbox_showoutoffov, handles.checkbox_shownonmatched, handles.checkbox_showmatched],'Units','pixels')



% Set Figure size and position
d=get(0,'screensize');
handles.uifigure_regrois.Position = [0.1*d(3), 0.1*d(4), 0.7*d(3), 0.7*d(4)];
% This automatically calls the changesize function


% Setup palette for plots (needs to be saved in handles)...
% Color is per session, linespec is per roi type - matched, nonmatched,
% out_of_fov/deleted are for now combined (see 'key2linespec')
handles.palette.tiff.colors = struct(...
    'unselected', {[50,150,240]/255, [230,157,0]/255}, ...
    'selected', {[232,13,110]/255, [20,200,100]/255});
handles.palette.tiff.linespec = struct(...
    'linewidth', {1,2,2}, ...
    'linewidth_select', {2,3,3}, ...
    'linestyle', {'-',':',':'});

% Tags definition ( matched, nonmatched, outoffov, deleted)
handles.keys.tags_def = {[0,1,2], [-1,-2], -3, -4};

% Initialise registration index to first one
handles.keys.current_registration = 1;
% Get the corresponding sessions indices in REG_DATA
[~, handles.keys.sessions_data_idx] = ismember(handles.REG_RES(handles.keys.current_registration).sessions, {handles.REG_DATA.session});
% Get the session allocation in REG_RES
handles.keys.sessions_reg_idx = [1,2];
% Get the name of the sessions
handles.keys.sesnames = handles.REG_RES(handles.keys.current_registration).sessions;
% Set default is 'aligned' to session 1
handles.keys.is_ref = [true,false];  

% Save session / axis allocation
handles.keys.sessions_uiaxes = {'uiaxes1', 'uiaxes2'};

% Show the sessions names for the radiobuttons to swap ref session
set(handles.radiobutton_s1, 'String', handles.keys.sesnames(1))
set(handles.radiobutton_s2, 'String', handles.keys.sesnames(2))

% Axes are not linked by default, so disable 'mask aligned to'
set([handles.radiobutton_s1, handles.radiobutton_s2], 'Enable', 'off')

% Add a convenient list of ROIs nums as strings...
[handles.keys.str_ROIsIDs{1,1}, handles.keys.n1n2(1,1)] = convert_ROIsIDs_to_strings(handles.REG_RES(handles.keys.current_registration), handles.keys.sessions_reg_idx(1));
[handles.keys.str_ROIsIDs{1,2}, handles.keys.n1n2(1,2)] = convert_ROIsIDs_to_strings(handles.REG_RES(handles.keys.current_registration), handles.keys.sessions_reg_idx(2));

% Initialise popupmenus Session
list_sessions = {handles.REG_DATA.session};
list_ses2 = setdiff(list_sessions, handles.keys.sesnames{1}, 'stable');
list_ses1 = setdiff(list_sessions, handles.keys.sesnames{2}, 'stable');
set(handles.popupmenu_s1, 'String', list_ses1, 'Value', find(strcmp(list_ses1,handles.keys.sesnames{1})))
set(handles.popupmenu_s2, 'String', list_ses2, 'Value', find(strcmp(list_ses2,handles.keys.sesnames{2})))

% Initialise popupmenus Inspect rois
list_inspect_type = {'Matched Pairs', 'Conflicting Chain', 'Incomplete Chain', 'Missed Pairs', 'Nonmatched 1', 'Nonmatched 2', 'Out of Fov 1', 'Out of Fov 2', 'Deleted rois 1', 'Deleted rois 2'};
set(handles.popupmenu_inspect, 'String', list_inspect_type)
% Save the list indices
handles.keys.inspect_type_idx = struct('matched_pairs', 1, 'conflicting_chain', 2, 'incomplete_chain', 3, 'missed_pairs', 4, 'nonmatched_ROIs', [5,6], 'outofFOV_ROIs', [7,8], 'deleted_ROIs', [9,10]);
handles.keys.inspect_type_name = {'matched_pairs', 'conflicting_chain', 'incomplete_chain', 'missed_pairs', 'nonmatched_ROIs', 'nonmatched_ROIs', 'outofFOV_ROIs', 'outofFOV_ROIs', 'deleted_ROIs', 'deleted_ROIs'};
handles.gui_state.popup_pairtype = [];

% Colors of list of rois
set(handles.listbox_rois_s1, 'BackgroundColor', [50,150,240]/255)
set(handles.listbox_rois_s2, 'BackgroundColor', [230,157,0]/255)

% Initialise SESSIONS: rois tags (matched, unmatched, deleted, etc...),
% graphs and List of rois
[handles] = initialise_sessions_gui(handles, true, 1, 'pointer', [1,1]);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_regROIs wait for user response (see UIRESUME)
% uiwait(handles.uifigure_regrois);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_regROIs_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = [];

% --- Executes when user attempts to close uifigure_regrois.
function uifigure_regrois_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to uifigure_regrois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

if ~isempty(handles)
    assignin('base', handles.output_var, handles.REG_RES)
end

delete(hObject);


% --- Executes when uifigure_regrois is resized.
function uifigure_regrois_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to uifigure_regrois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

the_position = handles.uifigure_regrois.Position;

% Buttons width and height (in pixels)
buttons_siz.w = 110;
buttons_siz.h = 25;
buttonscheck_siz.w = 120;
buttonscheck_siz.h = buttons_siz.h;
popupmenu_h = 22;

% Margins for objects with UI frame, in pixels
pix_margins.w = 5;
pix_margins.h = 2;
pix_margins.dw_2imgs = 5; % pixels spacing between images
pix_margins.dh_imgpa = 10; % pixels spacing between image and lower panel
% Margins for objects with frames of panels
pix_margins.wp = 5; % margin for width
pix_margins.hp_top = 20; % margin for content at top of panel
pix_margins.hp_bottom = 2; % margin for content at bottom
pix_margins.dh_buttons_min = 2; % minimum of 2pix between buttons
pix_margins.dh_buttons_max = 10; % max of 2pix between buttons

% Proportions
p_images=0.85; % total width for the images, proportional to total figure width

% Compute lower panel minimum size - based on Edit lower panel
h_min = pix_margins.hp_top + 3*buttons_siz.h + 2*pix_margins.dh_buttons_min + pix_margins.hp_bottom;

% Compute size of a tiff axis, taking into account the min panel size
image_wh_max = the_position(4) - h_min - 2*pix_margins.h - pix_margins.dh_imgpa;
image_wh = (the_position(3)*p_images - pix_margins.w - pix_margins.dw_2imgs)/2; % width for 1 tiff axis
image_wh = min(image_wh, image_wh_max);

% Resize Image axes
x_ax1 = pix_margins.w;
x_ax2 = x_ax1 + image_wh + pix_margins.dw_2imgs;
y = the_position(4) - pix_margins.h - image_wh;
handles.uiaxes1.Position = [x_ax1,y,image_wh,image_wh];
handles.uiaxes2.Position = [x_ax2,y,image_wh,image_wh];

% *** Upper pannel inspection ***
% Panel 
panel_h_upper = image_wh + pix_margins.h;
y = the_position(4) - panel_h_upper;
w = the_position(3)*(1-p_images) - pix_margins.w;
h = panel_h_upper - pix_margins.h;
x = the_position(3) - pix_margins.w - w;
handles.uipanel_inspect.Position = [x, y, w, h];
% Content -
dx = 1;
y = h-50;
handles.popupmenu_inspect.Position = [dx, y, w-dx, popupmenu_h];
dx_in = 3;
y = y-popupmenu_h-10;
w0 = w/2 - dx -dx_in/2;
handles.popupmenu_s1.Position = [dx, y, w0, popupmenu_h];
handles.popupmenu_s2.Position = [dx+w0+dx_in, y, w0, popupmenu_h];
y = 1;
h0 = handles.popupmenu_s1.Position(2) - 10 -1;
handles.listbox_rois_s1.Position = [dx, y, w0, h0];
handles.listbox_rois_s2.Position = [dx+w0+dx_in, y, w0, h0];

% *** Lower panels ***
p_panels = [0.4,0.25,0.35]; % proportions

% *** Lower panel View ***
panel_w = p_panels(1)*the_position(3) - pix_margins.w;
panel_h_lower = the_position(4) - panel_h_upper - pix_margins.dh_imgpa - pix_margins.h;
handles.uipanel_view.Position = [pix_margins.w, pix_margins.h, panel_w, panel_h_lower];
% Content - buttons
%   > reset view (top)
button_y = panel_h_lower - buttons_siz.h - pix_margins.hp_top;
handles.pushbutton_resetview.Position = [pix_margins.wp, button_y, buttons_siz.w, buttons_siz.h];
%   > compute spacing between buttons
dy_buttons = (panel_h_lower - pix_margins.hp_top - pix_margins.hp_bottom - 3*buttons_siz.h)/2;
dy_buttons = min(pix_margins.dh_buttons_max, dy_buttons);
%   > hide tiff 1 & 2
button_y = button_y - buttons_siz.h - dy_buttons;
handles.togglebutton_hidetiff1.Position = [pix_margins.wp, button_y, buttons_siz.w, buttons_siz.h];
button_y = button_y - buttons_siz.h - dy_buttons;
handles.togglebutton_hidetiff2.Position = [pix_margins.wp, button_y, buttons_siz.w, buttons_siz.h];
%   > Link axes
text_x = pix_margins.wp + handles.togglebutton_hidetiff2.Position(3);
text_y = panel_h_lower - pix_margins.hp_top - handles.text_linkaxes.Position(4);
handles.text_linkaxes.Position = [text_x, text_y, 60, handles.text_linkaxes.Position(4)];
button_y = text_y - handles.checkbox_linkaxes.Position(4) - 2;
button_x = text_x + 25;
handles.checkbox_linkaxes.Position = [button_x, button_y, handles.checkbox_linkaxes.Position(3:4)];
%   > subpanel select session ref
subp_h = panel_h_lower - pix_margins.hp_top;
subp_y = 0;
subp_x = handles.text_linkaxes.Position(1) + handles.text_linkaxes.Position(3);
subp_w = panel_w - subp_x - pix_margins.wp;
handles.uibuttongroup_alignmask.Position = [subp_x, subp_y, subp_w, subp_h];
button_x = 5;
button_w = subp_w-button_x-5;
button_y = subp_h-handles.radiobutton_s1.Position(4) - 30;
handles.radiobutton_s1.Position = [button_x, button_y, button_w, handles.radiobutton_s1.Position(4)];
button_y = button_y - handles.radiobutton_s2.Position(4) - 10;
handles.radiobutton_s2.Position = [button_x, button_y, button_w, handles.radiobutton_s2.Position(4)];


% *** Lower panel Show Rois ***
panel_x = handles.uipanel_view.Position(1) + panel_w;
panel_w = p_panels(2)*the_position(3);
handles.uipanel_showrois.Position = [panel_x, pix_margins.h, panel_w, panel_h_lower];

dx = (panel_w - 2*buttonscheck_siz.w)/3.5;
button_x = dx;
handles.checkbox_showmatched.Position = [button_x, handles.pushbutton_resetview.Position(2), buttonscheck_siz.w, buttonscheck_siz.h];
handles.checkbox_shownonmatched.Position = [button_x, handles.togglebutton_hidetiff1.Position(2), buttonscheck_siz.w, buttonscheck_siz.h];
handles.checkbox_showoutoffov.Position = [button_x, handles.togglebutton_hidetiff2.Position(2), buttonscheck_siz.w, buttonscheck_siz.h];

button_x = button_x + buttonscheck_siz.w + dx;
handles.checkbox_showdeleted.Position = [button_x, handles.checkbox_showmatched.Position(2), buttonscheck_siz.w, buttonscheck_siz.h];

text_x = button_x + handles.checkbox_showboth.Position(3)-2;
text_y = handles.checkbox_shownonmatched.Position(2) + buttonscheck_siz.h - handles.text_showboth.Position(4);
handles.text_showboth.Position = [text_x, text_y, handles.text_showboth.Position(3:4)];

button_y = handles.checkbox_shownonmatched.Position(2);
handles.checkbox_showboth.Position = [button_x, button_y, handles.checkbox_showboth.Position(3), buttonscheck_siz.h];



% *** Lower panel Edit Results ***
panel_x = handles.uipanel_showrois.Position(1) + panel_w;
panel_w = p_panels(3)*the_position(3) - pix_margins.w;
handles.uipanel_editresults.Position = [panel_x, pix_margins.h, panel_w, panel_h_lower];
% Content -
p_content = [0.3,0.7]; % proportions
box_h = 30;
text_x = (p_content(1)*panel_w - handles.text_distance.Position(3))/2;
text_y = handles.togglebutton_hidetiff2.Position(2) + (panel_h_lower - pix_margins.hp_top - handles.togglebutton_hidetiff2.Position(2))/2;
handles.text_distance.Position = [text_x, text_y, handles.text_distance.Position(3:4)];
box_w = handles.text_distance.Position(3); 
box_x = text_x;
box_y = text_y - box_h;
handles.text_distancevalue.Position = [box_x, box_y, box_w, box_h];

dx = (p_content(2)*panel_w - 2*buttons_siz.w)/3;
button_x = p_content(1)*panel_w + dx;
handles.pushbutton_match.Position = [button_x, handles.pushbutton_resetview.Position(2), buttons_siz.w, buttons_siz.h];
button_y = handles.togglebutton_hidetiff2.Position(2);
handles.pushbutton_deletepair.Position = [button_x, button_y, buttons_siz.w, buttons_siz.h];

button_x = handles.pushbutton_match.Position(1) + buttons_siz.w +dx;
handles.pushbutton_addroi.Position = [button_x, handles.pushbutton_match.Position(2), buttons_siz.w, buttons_siz.h];
handles.pushbutton_outoffov.Position = [button_x, handles.togglebutton_hidetiff1.Position(2), buttons_siz.w, buttons_siz.h];
handles.pushbutton_deleteroi.Position = [button_x, handles.togglebutton_hidetiff2.Position(2), buttons_siz.w, buttons_siz.h];







%__________________________________________________________________________
%
%   VARIOUS UTILITIES
%__________________________________________________________________________


function [r] = get_roi_list_from_pair_list(r_list, k)
if isempty(r_list)
    r=[];
else
    r = r_list(k,:);
end

function [str_ROIsIDs, n] = convert_ROIsIDs_to_strings(REG_RES, session_reg_idx)

% Add a convenient list of ROIs nums as strings...
[n] = length(REG_RES.ROIsIDs{session_reg_idx});
str_ROIsIDs = arrayfun(@(k)num2str(k), REG_RES.ROIsIDs{session_reg_idx}, 'UniformOutput', false);

function [img, imgoutline, rois] = get_session_data(handles, axes_are_linked)

img = cell(1,2);
imgoutline = cell(1,2);
rois = cell(1,2);


if axes_are_linked

    ref_ses = handles.keys.sessions_data_idx(handles.keys.is_ref);
    non_ref_ses = handles.keys.sessions_data_idx(~handles.keys.is_ref);
    
    
    % Get the ref and aligned sync tiff image
    img{handles.keys.is_ref} = handles.REG_DATA(ref_ses).aligned_data.(['s',handles.keys.sesnames{~handles.keys.is_ref}]).imgtiff.original_sync;
    img{~handles.keys.is_ref} = handles.REG_DATA(non_ref_ses).aligned_data.(['s',handles.keys.sesnames{handles.keys.is_ref}]).imgtiff.aligned_sync;
    
    % Prep ROIsXY
    u=cell(1,2);
    u{handles.keys.is_ref} = 'original_sync';
    
    u{~handles.keys.is_ref} = 'aligned_sync';
    rois{handles.keys.is_ref} = { handles.REG_DATA(handles.keys.sessions_data_idx(1)).aligned_data.(['s',handles.keys.sesnames{2}]).ROIs_XY.(u{1}), ...
        handles.REG_DATA(handles.keys.sessions_data_idx(2)).aligned_data.(['s',handles.keys.sesnames{1}]).ROIs_XY.(u{2}) };
    
    u{~handles.keys.is_ref} = 'globalT_sync';
    rois{~handles.keys.is_ref} = { handles.REG_DATA(handles.keys.sessions_data_idx(1)).aligned_data.(['s',handles.keys.sesnames{2}]).ROIs_XY.(u{1}), ...
        handles.REG_DATA(handles.keys.sessions_data_idx(2)).aligned_data.(['s',handles.keys.sesnames{1}]).ROIs_XY.(u{2}) };
    
    % Prep image outline
    u=cell(1,2);
    u{handles.keys.is_ref} = 'original_sync';
    u{~handles.keys.is_ref} = 'aligned_sync';
    
    imgoutline{handles.keys.is_ref} = { handles.REG_DATA(handles.keys.sessions_data_idx(1)).aligned_data.(['s',handles.keys.sesnames{2}]).corners_XY.(u{1}), ...
        handles.REG_DATA(handles.keys.sessions_data_idx(2)).aligned_data.(['s',handles.keys.sesnames{1}]).corners_XY.(u{2}) };
    
    imgoutline{~handles.keys.is_ref} = { handles.REG_DATA(handles.keys.sessions_data_idx(1)).aligned_data.(['s',handles.keys.sesnames{2}]).corners_XY.(u{1}), ...
        handles.REG_DATA(handles.keys.sessions_data_idx(2)).aligned_data.(['s',handles.keys.sesnames{1}]).corners_XY.(u{2}) };
    

else
    
    % Keep template images
    img{1} = handles.REG_DATA(handles.keys.sessions_data_idx(1)).data.imgtiff;
    img{2} = handles.REG_DATA(handles.keys.sessions_data_idx(2)).data.imgtiff;
    
    % Prep ROIsXY
    rois{1} = { handles.REG_DATA(handles.keys.sessions_data_idx(1)).data.ROIs_XY , ...
        handles.REG_DATA(handles.keys.sessions_data_idx(2)).aligned_data.(['s',handles.keys.sesnames{1}]).ROIs_XY.aligned };
    rois{2} = { handles.REG_DATA(handles.keys.sessions_data_idx(1)).aligned_data.(['s',handles.keys.sesnames{2}]).ROIs_XY.aligned, ...
        handles.REG_DATA(handles.keys.sessions_data_idx(2)).data.ROIs_XY };
    
    % Prep image outline
    
    imgoutline{1} = { handles.REG_DATA(handles.keys.sessions_data_idx(1)).data.corners_XY, ...
        handles.REG_DATA(handles.keys.sessions_data_idx(2)).aligned_data.(['s',handles.keys.sesnames{1}]).corners_XY.aligned };
    
    imgoutline{2} = { handles.REG_DATA(handles.keys.sessions_data_idx(1)).aligned_data.(['s',handles.keys.sesnames{2}]).corners_XY.aligned, ...
        handles.REG_DATA(handles.keys.sessions_data_idx(2)).data.corners_XY };
       
    
end


%__________________________________________________________________________
%
%   UTILITIES FOR ROIs LISTS
%__________________________________________________________________________


function [min_c_dist] = get_closest_centroids_distance(c_ref, c_test)

% Get the min centroids distances for each ref
min_c_dist = cellfun(@(c0) min(cellfun(@(c) sum((c0-c).^2),c_test)), c_ref);

function [idx_c_select] = pick_closest_centroids(c_ref, c_test)

n_select = min(length(c_test), 5);

% Get the centroids distances
c_dist = cellfun(@(c) sum((c_ref-c).^2), c_test);

% Get the 5 closest centroids
[~, idx_c] = sort(c_dist, 'ascend');
idx_c_select = idx_c(1:n_select);

function [rois_nonmatched] = pick_closest_nonmatched_rois(list_nonmatched2update, c_ref, c_test)

if isempty(list_nonmatched2update)
    rois_nonmatched = [];
else
    % Get the ROIs with closest centroids
    [idx_c_select] = pick_closest_centroids(c_ref, c_test);
    rois_nonmatched = list_nonmatched2update(idx_c_select);
    
end

function [selected_roi] = update_rois_list(h_listROIs, rois, str_ROIsIDs, pointer)


% Check if allow selection
if isempty(rois)
    set(h_listROIs, 'Enable', 'off', 'String', [])
    selected_roi = [];
else
    % Show the new list
    set(h_listROIs, 'Enable', 'on', 'String', str_ROIsIDs(rois))
    
    pointer = min(pointer, length(rois));
    
    set(h_listROIs, 'Value', pointer)
    selected_roi = rois(pointer);
end


%__________________________________________________________________________
%
%   PLOT UTILITIES
%__________________________________________________________________________

function [ls] = roistag2linespec(rois_tag)
ls = 2*ones(size(rois_tag)); ls(rois_tag>=0)=1; ls(rois_tag<-2)=3;

function plot_ROI(the_UIAxe_handle, ROIs_XY, c, lw, ls)

plot(the_UIAxe_handle, [ROIs_XY(:,1); ROIs_XY(1,1)], [ROIs_XY(:,2); ROIs_XY(1,2)], ...
    'Color', c, 'LineWidth', lw, 'LineStyle', ls)

function draw_image_outline(the_UIAxe_handle, C_XY, the_color)
% Note: corners_XY_zzz : [top left, bottom left, bottom right, top right]

plot(the_UIAxe_handle, [C_XY(1,:), C_XY(1,1)],[C_XY(2,:), C_XY(2,1)],'LineWidth',1,'Color', the_color)

function change_ROIs_masks(h_child_mask,ROIs_XY)

% Shift XY coordinates rois masks data 
for k=1:length(ROIs_XY)
    
    h_child_mask(k).XData = [ ROIs_XY{k}(:,1)', ROIs_XY{k}(1,1)];
    
    h_child_mask(k).YData = [ ROIs_XY{k}(:,2)', ROIs_XY{k}(1,2)];

end

function change_image_outline(h_child_imgo,C_XY)

% Shift XY coordinates of the frame
if iscell(C_XY)
    for k=1:2
        h_child_imgo(k).XData = [C_XY{k}(1,:), C_XY{k}(1,1)];
        h_child_imgo(k).YData = [C_XY{k}(2,:), C_XY{k}(2,1)];
    end
else
    h_child_imgo.XData = [C_XY(1,:), C_XY(1,1)];
    h_child_imgo.YData = [C_XY(2,:), C_XY(2,1)];
end

function update_ROIs_highlight(h_child_mask, previous_roi, new_roi, rois_tag, pal_color, linespec)

if ~isempty(previous_roi)
    ls_idx = roistag2linespec(rois_tag(previous_roi));
    % Deselect
    h_child_mask(previous_roi).LineWidth = linespec(ls_idx).linewidth;
    h_child_mask(previous_roi).Color = pal_color.unselected; 
end

if ~isempty(new_roi)
    ls_idx = roistag2linespec(rois_tag(new_roi));
    % Select
    h_child_mask(new_roi).LineWidth = linespec(ls_idx).linewidth_select;
    h_child_mask(new_roi).Color = pal_color.selected;
end

function swap_single_mask_roi_type(h_axes, roi, linespec)

for k=1:length(h_axes)
    h_axes(k).Children(roi).LineStyle = linespec.linestyle;
end

function swap_multi_masks_rois_type(h_axes, rois, rois_tag, the_linespec)

ls_idx = roistag2linespec(rois_tag);

for k=1:length(h_axes)
    
    for iRoi=1:length(rois)
        h_axes(k).Children(rois(iRoi)).LineStyle = the_linespec(ls_idx(iRoi)).linestyle;
    end
    
end


%__________________________________________________________________________
%
%   INITIALISATION FUNCTIONS
%__________________________________________________________________________

function [handles] = initialise_sessions_gui(handles, is_first_init, popup_pairtype, pointer_type, pointers)

% Adjust the value of pair inspection
set(handles.popupmenu_inspect, 'Value', popup_pairtype)

% Re-Initialise rois tags (matched, unmatched, deleted, etc...)
handles.keys.rois_tag = init_ROIs_tags(handles.REG_RES(handles.keys.current_registration), handles.keys.sessions_reg_idx);

% Get the data for the graphs
axes_are_linked = get(handles.checkbox_linkaxes,'Value')==1;
[img, imgoutline, rois] = get_session_data(handles, axes_are_linked);

% Initialise graphs
handles.orig_siz.plotlimits=cell(1,2);
handles.orig_siz.plotlimits{1} = init_uiaxe(handles, is_first_init, [true,false], img{1}, imgoutline{1}, rois{1});
handles.orig_siz.plotlimits{2} = init_uiaxe(handles, is_first_init, [false,true], img{2}, imgoutline{2}, rois{2});

if ~is_first_init
   checkbox_showboth_Callback(handles.checkbox_showboth, [], handles) 
   checkbox_showmatched_Callback(handles.checkbox_showmatched, [], handles) 
   checkbox_shownonmatched_Callback(handles.checkbox_shownonmatched, [], handles) 
   checkbox_showoutoffov_Callback(handles.checkbox_showoutoffov, [], handles) 
   checkbox_showdeleted_Callback(handles.checkbox_showdeleted, [], handles) 
end

% Get the largest limits
handles.orig_siz.plotlimits_max.Xlim = [min([handles.orig_siz.plotlimits{1}.Xlim(1), handles.orig_siz.plotlimits{2}.Xlim(1)]), ...
    max([handles.orig_siz.plotlimits{1}.Xlim(2), handles.orig_siz.plotlimits{2}.Xlim(2)])];
handles.orig_siz.plotlimits_max.Ylim = [min([handles.orig_siz.plotlimits{1}.Ylim(1), handles.orig_siz.plotlimits{2}.Ylim(1)]), ...
    max([handles.orig_siz.plotlimits{1}.Ylim(2), handles.orig_siz.plotlimits{2}.Ylim(2)])];

% Initialise List of rois
handles.keys.rois_pool = cell(1,2); handles.keys.selected_rois = cell(1,2);
[handles.keys.rois_pool, handles.keys.selected_rois] = get_rois_data_for_inspection(handles, [], popup_pairtype, pointer_type, pointers);
% Update current popup pair type choice
handles.gui_state.popup_pairtype = popup_pairtype;


function [plotlimits_init] = init_uiaxe(handles, is_first_init, is_active_ses, img, imgoutline, rois)

the_UIAxe = handles.keys.sessions_uiaxes{is_active_ses};

if is_first_init
    
    % Show tiff image
    axes(handles.(the_UIAxe))
    imagesc(img, [0,255]);
    colormap(handles.(the_UIAxe),gray);
    hold(handles.(the_UIAxe),'on')

    % Draw image outline
    draw_image_outline(handles.(the_UIAxe), imgoutline{2}, handles.palette.tiff.colors(2).unselected)
    draw_image_outline(handles.(the_UIAxe), imgoutline{1}, handles.palette.tiff.colors(1).unselected)
    % Hide the active session 
    h = handles.(the_UIAxe).Children(1:2);
    set(h(is_active_ses),'Visible', 'off')
    
else
    % Change the tiff images
    handles.(the_UIAxe).Children(end).CData = img;

    % Change the outlines image for the other session
    change_image_outline(handles.(the_UIAxe).Children(end-1),imgoutline{2})
    change_image_outline(handles.(the_UIAxe).Children(end-2),imgoutline{1})
    
    % Delete the ROIs
    delete(handles.(the_UIAxe).Children(1:end-3)) % at the end we have: image outline x2 and tiff image   
    
end

% Make boundaries tight to image
axis(handles.(the_UIAxe), 'tight')

% Draw the ROIs

% First session 2
init_ROIs_masks(handles.(the_UIAxe), rois{2}, handles.keys.rois_tag{2}, handles.palette.tiff.colors(2).unselected, handles.palette.tiff.linespec)
% Then session 1
init_ROIs_masks(handles.(the_UIAxe), rois{1}, handles.keys.rois_tag{1}, handles.palette.tiff.colors(1).unselected, handles.palette.tiff.linespec)

% Finalise axes
set(handles.(the_UIAxe), 'XTick',[], 'YTick', [], 'Color', [0,0,0])
axis(handles.(the_UIAxe), 'equal')
% Save limits
plotlimits_init.Xlim = get(handles.(the_UIAxe), 'Xlim');
plotlimits_init.Ylim = get(handles.(the_UIAxe), 'Ylim');


        
function init_ROIs_masks(h_UIAxe, ROIs_XY, rois_tag, the_color, the_linespec)
% Invert the handles of rois so that matches 1->Nrois in the handles.
% (remember: last is image, then have the image outline of other ses, and
% then stacks the rois masks of both sessions ; so the first roi is at n-2
% for session 2, then session 1.

ls_idx = roistag2linespec(rois_tag);

for k=length(ROIs_XY):-1:1 
        
    plot_ROI(h_UIAxe, ROIs_XY{k}, the_color, the_linespec(ls_idx(k)).linewidth, the_linespec(ls_idx(k)).linestyle)
        
end
    

function [rois_tags] = init_ROIs_tags(REG_RES, sessions_reg_idx)

rois_tags = cell(1,2);

for reg_iSes=1:2
    
    list_all = [...
        get_roi_list_from_pair_list(REG_RES.matched_pairs, reg_iSes), ...
        get_roi_list_from_pair_list(REG_RES.conflicting_chain, reg_iSes), ...
        get_roi_list_from_pair_list(REG_RES.incomplete_chain, reg_iSes), ...
        REG_RES.nonmatched_ROIs{reg_iSes}, ...
        get_roi_list_from_pair_list(REG_RES.missed_pairs, reg_iSes), ...
        REG_RES.outofFOV_ROIs{reg_iSes}, REG_RES.deleted_ROIs{reg_iSes}];
    
    [~,idx_sort] = sort(list_all, 'ascend');
    
    rois_tags{reg_iSes} = [...
        zeros(1,size(REG_RES.matched_pairs,2)), ...
        ones(1,size(REG_RES.conflicting_chain,2)), ...
        2*ones(1, size(REG_RES.incomplete_chain, 2)), ...
        -ones(1,length(REG_RES.nonmatched_ROIs{reg_iSes})), ...
        -2*ones(1,size(REG_RES.missed_pairs,2)), ...
        -3*ones(1,length(REG_RES.outofFOV_ROIs{reg_iSes})), ...
        -4*ones(1,length(REG_RES.deleted_ROIs{reg_iSes}))];
    
    % Match tags to original ROIs order
    rois_tags{reg_iSes} = rois_tags{reg_iSes}(idx_sort);
    
end
rois_tags = rois_tags(sessions_reg_idx);

   

%__________________________________________________________________________
% 
%        POPUPMENU SESSIONS
%           - SWAP REGISTRATION INSPECTED ON GUI
%__________________________________________________________________________

% --- Executes on selection change in popupmenu_s1.
function popupmenu_s1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_s1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_s1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_s1

contents = cellstr(get(hObject,'String'));
new_sessions = contents{get(hObject,'Value')};
% First check actually selected a new session
if ~strcmp(new_sessions, handles.keys.sesnames{1})
    swap_registration(handles, [true, false], new_sessions)
end

% --- Executes on selection change in popupmenu_s2.
function popupmenu_s2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_s2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_s2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_s2

contents = cellstr(get(hObject,'String'));
new_sessions = contents{get(hObject,'Value')};
% First check actually selected a new session
if ~strcmp(new_sessions, handles.keys.sesnames{2})
    swap_registration(handles, [false, true], new_sessions)
end

function swap_registration(handles, is_updated_ses, new_session)
  
% Get the corresponding sessions indices in REG_DATA
new_sessions_data_idx = find(strcmp(new_session, {handles.REG_DATA.session}));
% New pair of sesnames
new_sessions = handles.keys.sesnames;
new_sessions{is_updated_ses} = new_session;

all_paired_session = {handles.REG_RES.sessions};
% Find the new registration index corresponding to the new pair of sessions
new_registration = find(cellfun(@(x)all(ismember(new_sessions,x)),all_paired_session,'un', true));

% Get the corresponding sessions indices in REG_RES of current registration
[~, new_keys_sessions_reg_idx] = ismember(new_sessions, handles.REG_RES(new_registration).sessions);

% Save new registration index
handles.keys.current_registration = new_registration;
% Save the name of the new session
handles.keys.sesnames = new_sessions;
% Get the corresponding sessions indices in REG_DATA
handles.keys.sessions_data_idx(is_updated_ses) = new_sessions_data_idx;
% Get the session allocation in REG_RES
handles.keys.sessions_reg_idx = new_keys_sessions_reg_idx;

% Change the session name for the radiobuttons to swap ref session
h_radiobutton_sX = [handles.radiobutton_s1, handles.radiobutton_s2];
set(h_radiobutton_sX(is_updated_ses), 'String', new_session)

% Re-Initialise popupmenu for the session that was not updated
h_popupmenu_sX = [handles.popupmenu_s1, handles.popupmenu_s2];
list_sessions = {handles.REG_DATA.session};
list_ses = setdiff(list_sessions, new_session, 'stable');
set(h_popupmenu_sX(~is_updated_ses), 'String', list_ses, 'Value', find(strcmp(list_ses,new_sessions{~is_updated_ses})))

% Change the convenient list of ROIs nums as strings...
[handles.keys.str_ROIsIDs{is_updated_ses}, handles.keys.n1n2(is_updated_ses)] = convert_ROIsIDs_to_strings(handles.REG_RES(new_registration), new_keys_sessions_reg_idx(is_updated_ses));

% Track the roi of the session that was not updated
tracked_roi = handles.keys.selected_rois{~is_updated_ses};

if isempty(tracked_roi)
    popup_pairtype = handles.gui_state.popup_pairtype;
    pointers = [1,1];
    pointer_type = 'pointer';
else
    % Find it in the new registration
    S = fieldnames(handles.keys.inspect_type_idx);
    roi_not_found = true; k=0;
    pointers = [1,1];
    pointer_type = 'roi_tracker';
    while roi_not_found
        k=k+1;
        switch S{k}
            case {'matched_pairs','conflicting_chain','incomplete_chain','missed_pairs'}
                if ~isempty(handles.REG_RES(new_registration).(S{k}))
                    has_roi = tracked_roi == handles.REG_RES(new_registration).(S{k})(new_keys_sessions_reg_idx(~is_updated_ses),:);
                    roi_not_found = ~any(has_roi);
                    if ~roi_not_found
                        % Get the new popup menu for type inspection
                        popup_pairtype = k;
                        pointers(~is_updated_ses) = tracked_roi;
                        pointers(is_updated_ses) = handles.REG_RES(new_registration).(S{k})(new_keys_sessions_reg_idx(is_updated_ses),has_roi);
                    end
                end
            otherwise
                if ~isempty(handles.REG_RES(new_registration).(S{k}){new_keys_sessions_reg_idx(~is_updated_ses)})
                    roi_not_found = ~any(tracked_roi == handles.REG_RES(new_registration).(S{k}){new_keys_sessions_reg_idx(~is_updated_ses)});
                    if ~roi_not_found
                        % Get the new popup menu for type inspection
                        popup_pairtype = handles.keys.inspect_type_idx.(S{k})(~is_updated_ses);
                        pointers(~is_updated_ses) = tracked_roi;
                    end
                end
        end
    end
end

[handles] = initialise_sessions_gui(handles, false, popup_pairtype, pointer_type, pointers);

% Update handles structure
guidata(handles.uifigure_regrois, handles);


%__________________________________________________________________________
% 
%        POPUPMENU & ROIs LISTS
%           - CHOOSE TYPE TO INSPECT: PAIRED, NON-MATCHED etc.
%           - INSPECT ROIS... INTERACT WITH LIST
%           -> highligh selected rois and show pair distance
%__________________________________________________________________________
        

function [rois, rois_selected] = get_rois_data_for_inspection(handles, previous_popup_pairtype, popup_pairtype, varargin)

% varargin : pointer to selected rois 
% -> If it is an internal call (from delete/pair buttons), in that case
% keep the pointer where it is so that automatically select the next item
% in the list. 
% -> If change session, we track the roi of the session that was not
% updated, and therefore change the popupmenu for inspection type. For any
% non-paired elements, the second pointer has to be set to 1 by default.
% -> Otherwise reset the pointer at the start of the list.
if isempty(varargin)
    pointers = [1,1];
    check_pointers_track_rois = false;
else
    pointers = varargin{2};
    check_pointers_track_rois = strcmp(varargin{1}, 'roi_tracker');
end
    
rois = [];
rois_selected = [];
% Check that the user requested a different one! Otherwise do nothing
if isempty(previous_popup_pairtype) || previous_popup_pairtype~=popup_pairtype
    
    rois = cell(1,2);
    rois_selected = cell(1,2);

    inspect_type_name = handles.keys.inspect_type_name{popup_pairtype};
    
    if any( popup_pairtype == ...
            [handles.keys.inspect_type_idx.matched_pairs, handles.keys.inspect_type_idx.conflicting_chain, ...
            handles.keys.inspect_type_idx.incomplete_chain, handles.keys.inspect_type_idx.missed_pairs] ) % ALL POSSIBLE MATCHED ROIS
        
        if isempty(handles.REG_RES(handles.keys.current_registration).(inspect_type_name))
            rois{1} = []; rois{2} = [];
        else
            [n1,n2] = size(handles.REG_RES(handles.keys.current_registration).score);
            % Get the indices
            rois{1} = handles.REG_RES(handles.keys.current_registration).(inspect_type_name)(1, :);
            rois{2} = handles.REG_RES(handles.keys.current_registration).(inspect_type_name)(2, :);
            idx_pair = sub2ind([n1,n2], rois{1}, rois{2});
            % Sort by score (lower the better == distance)
            pair_distance = handles.REG_RES(handles.keys.current_registration).score(idx_pair);
            [~, idx_sort] = sort(pair_distance, 'descend');
            rois{1} = rois{1}(idx_sort);
            rois{2} = rois{2}(idx_sort);
            
            % We have to sort according to which session is 1 and 2
            rois = rois(handles.keys.sessions_reg_idx);
            
            % Set the pointers
            if check_pointers_track_rois
                pointers = [find(rois{1}==pointers(1)), find(rois{2}==pointers(2))];                
            end
        end
        
    elseif any( popup_pairtype == handles.keys.inspect_type_idx.nonmatched_ROIs ) % NON-MATCHED ROIS
        
        is_caller = popup_pairtype == handles.keys.inspect_type_idx.nonmatched_ROIs ;
        
        is_caller_ses_reg = is_caller(handles.keys.sessions_reg_idx);
        caller_ses_data_idx = handles.keys.sessions_data_idx(is_caller);
        noncaller_ses_data_idx = handles.keys.sessions_data_idx(~is_caller);
        
        if isempty(handles.REG_RES(handles.keys.current_registration).nonmatched_ROIs{is_caller})
            rois{1} = [];
            rois{2} = [];
        else
            % Get the rois indices of the inspected caller
            rois{is_caller} = handles.REG_RES(handles.keys.current_registration).nonmatched_ROIs{is_caller_ses_reg};
            rois_noncaller = handles.REG_RES(handles.keys.current_registration).nonmatched_ROIs{~is_caller_ses_reg};
            
            img_type = {'original','aligned'};
            
            % Sort by centroid distance so first display nonmatched rois
            % that are very close to other nonmatched rois
            caller_sesname = ['s',handles.keys.sesnames{is_caller}];
            noncaller_sesname = ['s',handles.keys.sesnames{~is_caller}];
            [min_c_dist] = get_closest_centroids_distance(handles.REG_DATA(caller_ses_data_idx).aligned_data.(noncaller_sesname).centroids.(img_type{is_caller})(rois{is_caller}), ...
                handles.REG_DATA(noncaller_ses_data_idx).aligned_data.(caller_sesname).centroids.(img_type{~is_caller})(rois_noncaller));
            [~, idx_sort] = sort(min_c_dist, 'ascend');
            rois{is_caller} = rois{is_caller}(idx_sort);
            
            % Set the pointers
            if check_pointers_track_rois
                pointers(is_caller) = find(rois{is_caller}==pointers(is_caller));
            end
            
            % Then select the top one automatically and Get the
            % corresponding list of rois from the other session
            [rois{~is_caller}] = pick_closest_nonmatched_rois(...
                handles.REG_RES(handles.keys.current_registration).nonmatched_ROIs{~is_caller_ses_reg}, ...
                handles.REG_DATA(caller_ses_data_idx).aligned_data.(noncaller_sesname).centroids.(img_type{is_caller}){rois{is_caller}(min(pointers(is_caller), length(rois{is_caller})))}, ...
                handles.REG_DATA(noncaller_ses_data_idx).aligned_data.(caller_sesname).centroids.(img_type{~is_caller})(rois_noncaller));
            
        end
        
    else % EDGES OR OUT OF FOVS ROIS / DELETED ROIS
        
        is_caller = popup_pairtype == handles.keys.inspect_type_idx.(inspect_type_name) ;
        is_caller_ses_reg = is_caller(handles.keys.sessions_reg_idx);
        
        % Get the indices
        rois{~is_caller} = [];
        rois{is_caller} = handles.REG_RES(handles.keys.current_registration).(inspect_type_name){is_caller_ses_reg};
        
        % Set the pointers
        if check_pointers_track_rois
            pointers(is_caller) = find(rois{is_caller}==pointers(is_caller));
        end
                
    end
    
    % Update lists
    rois_selected{1} = update_rois_list(handles.listbox_rois_s1, rois{1}, handles.keys.str_ROIsIDs{1}, pointers(1));
    rois_selected{2} = update_rois_list(handles.listbox_rois_s2, rois{2}, handles.keys.str_ROIsIDs{2}, pointers(2));
    

    % Update push buttons
    update_enable_push_buttons(popup_pairtype, cellfun(@(s)~isempty(s), rois_selected), handles.keys.inspect_type_idx, ...
        handles.pushbutton_deletepair, handles.pushbutton_match, handles.pushbutton_outoffov, handles.pushbutton_addroi, handles.pushbutton_deleteroi);
    
    % Update the distances / rois highlight on tiff axes
    update_ROISpair_on_gui(handles, [true, true], rois_selected)
    
    

end

% --- Executes on selection change in popupmenu_inspect.
function popupmenu_inspect_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_inspect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_inspect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_inspect

popup_pairtype = get(hObject, 'Value');

[rois, rois_selected] = get_rois_data_for_inspection(handles, handles.gui_state.popup_pairtype, popup_pairtype);

% Update current popup pair type choice
handles.gui_state.popup_pairtype = popup_pairtype;
    
% Check that the user requested a different one! Otherwise do nothing
if ~isempty(rois)
    % Save the rois keys
    handles.keys.selected_rois = rois_selected;
    handles.keys.rois_pool = rois;

    % Update handles structure
    guidata(handles.uifigure_regrois, handles);
end

% --- Executes on selection change in listbox_rois_s1.
function listbox_rois_s1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_rois_s1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_rois_s1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_rois_s1

interact_with_roi_list(handles,[true,false])

% --- Executes on selection change in listbox_rois_s2.
function listbox_rois_s2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_rois_s2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_rois_s2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_rois_s2

interact_with_roi_list(handles,[false,true])

function interact_with_roi_list(handles, is_caller)

h_listRois = {handles.listbox_rois_s1, handles.listbox_rois_s2};

% Get the roi calling...
new_value = get(h_listRois{is_caller}, 'Value');
roi_caller = handles.keys.rois_pool{is_caller}(new_value);
% Check if caller is unchaged...
if roi_caller == handles.keys.selected_rois{is_caller}
    % If so, do nothing
    new_rois = handles.keys.selected_rois;
    is_updated = false(1,2);
    new_rois_list = [];
else
    % Otherzise, operations depend on the the type of rois inspected...
    
    if any( handles.gui_state.popup_pairtype == ...
            [handles.keys.inspect_type_idx.matched_pairs, handles.keys.inspect_type_idx.conflicting_chain, ...
            handles.keys.inspect_type_idx.incomplete_chain, handles.keys.inspect_type_idx.missed_pairs] ) % ALL POSSIBLE MATCHED ROIS
        
        % ROIs are paired, hence change ROI POINTER for both sessions
        
        new_rois = cell(1,2);
        is_updated = true(1,2);
        new_rois_list = [];
        
        % Get the new roi
        new_rois{is_caller} = roi_caller;
        
        % Change the other session
        set(h_listRois{~is_caller}, 'Value', new_value)
        new_rois{~is_caller} = handles.keys.rois_pool{~is_caller}(new_value);
        
    elseif any( handles.gui_state.popup_pairtype == handles.keys.inspect_type_idx.nonmatched_ROIs ) % NON-MATCHED ROIS
        
        % ROIs are not paired, so
        % - if caller == nonmatch_ref , change ROIs for caller AND update
        % other session + reinit roi selection
        % - if caller ~= nonmatch_ref , change ROIs for caller only (i.e.
        % other session)
        
        % Get the session of the inspected roi
        is_ref = handles.gui_state.popup_pairtype == handles.keys.inspect_type_idx.nonmatched_ROIs ;
        
        % Check if caller == nonmatch_ref
        caller_is_ref = is_ref(is_caller);

        % Check if there are possible matches in the other session
        is_caller_ses_reg = is_caller(handles.keys.sessions_reg_idx);
        list_nonmatched2update = handles.REG_RES(handles.keys.current_registration).nonmatched_ROIs{~is_caller_ses_reg};
        has_possible_pairs = ~isempty(list_nonmatched2update);
        
        new_rois = cell(1,2);
        new_rois{is_caller} = roi_caller;
        
        if caller_is_ref && has_possible_pairs
            
            is_updated = true(1,2);
            
            % Update the new nonmatch of the other session, selecting the
            % closest ones based on centroids distance
            centroids = {handles.REG_DATA(handles.keys.sessions_data_idx(1)).aligned_data.(['s',handles.keys.sesnames{2}]).centroids.original, ...
                handles.REG_DATA(handles.keys.sessions_data_idx(2)).aligned_data.(['s',handles.keys.sesnames{1}]).centroids.aligned};
            str_ROIsIDs_2update = handles.keys.str_ROIsIDs{~is_caller};
            
            [new_rois_list] = pick_closest_nonmatched_rois(...
                list_nonmatched2update, ...
                centroids{is_caller}{roi_caller}, ...
                centroids{~is_caller}(list_nonmatched2update));
            
            % Update lists
            new_rois{~is_caller} = update_rois_list(h_listRois{~is_caller}, new_rois_list, str_ROIsIDs_2update, 1);
            
        else
            % If it is not the inspected roi, or if no possible match,
            % either way we update only one list: the caller
            is_updated = false(1,2);
            is_updated(is_caller) = true;
            new_rois(~is_caller) = handles.keys.selected_rois(~is_caller);
            new_rois_list = [];
            
        end
        
        
    else % EDGES OR OUT OF FOV, DELETED ROIS
        
        % Only update the caller
        new_rois = cell(1,2);
        is_updated = false(1,2);
        is_updated(is_caller) = true;
        new_rois_list = [];
        % Save the new roi
        new_rois{is_caller} = roi_caller;
    end
    
end

if any(is_updated)
    % Update the distances / rois highlight on tiff axes
    update_ROISpair_on_gui(handles, is_updated, new_rois)
    
    handles.keys.selected_rois = new_rois;
    
    if ~isempty(new_rois_list)
        handles.keys.rois_pool{~is_caller} = new_rois_list;
    end
    
    % Update handles structure
    guidata(handles.uifigure_regrois, handles);
end

function update_ROISpair_on_gui(handles, is_updated, new_rois)

% Update the distances / rois highlight on tiff axes

% We have to sort back who is session 1/2 in the GUI and in the
% registration data
new_rois_for_distance = [new_rois(handles.keys.sessions_reg_idx==1) , new_rois(handles.keys.sessions_reg_idx==2)];
% Update distance on GUI
update_pair_distance(handles.text_distancevalue, handles.REG_RES(handles.keys.current_registration).score, new_rois_for_distance)

% Update masks
idx_ses = 1;
if is_updated(idx_ses)
    
    % Tiffs
    update_ROIs_highlight(handles.uiaxes1.Children(1:handles.keys.n1n2(1)), ...
        handles.keys.selected_rois{idx_ses}, new_rois{idx_ses}, handles.keys.rois_tag{idx_ses}, ...
        handles.palette.tiff.colors(idx_ses), handles.palette.tiff.linespec)
    
    update_ROIs_highlight(handles.uiaxes2.Children(1:handles.keys.n1n2(1)), ...
        handles.keys.selected_rois{idx_ses}, new_rois{idx_ses}, handles.keys.rois_tag{idx_ses}, ...
        handles.palette.tiff.colors(idx_ses), handles.palette.tiff.linespec)
    
end

idx_ses = 2;
if is_updated(idx_ses)
    
    % Tiffs
    update_ROIs_highlight(handles.uiaxes2.Children(handles.keys.n1n2(1)+1:end-2), ...
        handles.keys.selected_rois{idx_ses}, new_rois{idx_ses}, handles.keys.rois_tag{idx_ses}, ...
        handles.palette.tiff.colors(idx_ses), handles.palette.tiff.linespec)
    update_ROIs_highlight(handles.uiaxes1.Children(handles.keys.n1n2(1)+1:end-2), ...
        handles.keys.selected_rois{idx_ses}, new_rois{idx_ses}, handles.keys.rois_tag{idx_ses}, ...
        handles.palette.tiff.colors(idx_ses), handles.palette.tiff.linespec)
    
end

function update_pair_distance(h_pair_distance, d, rois)

if isempty(rois{1}) || isempty(rois{2})
    set(h_pair_distance, 'String', [])
else
    set(h_pair_distance, 'String', num2str(round( d(rois{1}, rois{2}) *100)/100))
end

function update_enable_push_buttons(popup_pairtype, roi_is_selected, keys_inspect_type_idx, h_deletepair, h_match, h_outoffovroi, h_addroi, h_deleteroi)
 % EDGES OR OUT OF FOVS ROIS % DELETED ROIS

% Check current state :
% - if no selection, all disabled
% - if at least one selected, then check:
%       * if inspect pairs -> allow delete pair
%       * if inspect pairs from conflicting chains -> allow delete pair
%       * if inspect pairs from incomplete chains -> allow delete pair
%       * if inspect pairs from missing pairs -> allow match and out_of_fov / delete roi
%       * if inspect nonmatched -> allow match, out_of_fov and delete roi (current ses)
%       * if inspect out_of_fov -> allow add and delete roi (current ses)
%       * if inspect deleted -> allow add and out_of_fov roi (current ses)


        
list_h = {h_deletepair, h_match, h_outoffovroi, h_addroi, h_deleteroi};

% Default to turn off is if it is ON
is_OFF = cellfun(@(h) strcmp(get(h, 'Enable'), 'off'), list_h);
set_OFF = ~is_OFF;
set_ON = false(1,5);

if any(roi_is_selected)
    
    if any( popup_pairtype == [keys_inspect_type_idx.matched_pairs, keys_inspect_type_idx.conflicting_chain, ...
            keys_inspect_type_idx.incomplete_chain] ) % ALL MATCHED ROIS
        
        set_ON(1) = is_OFF(1); set_OFF(1) = false; % allow delete pair
        
    elseif popup_pairtype == keys_inspect_type_idx.missed_pairs % MISSING PAIRS
        
        set_ON(2) = is_OFF(2); set_OFF(2) = false; % allow create pair
        set_ON(3) = is_OFF(3); set_OFF(3) = false; % allow out_of_fov roi
        set_ON(5) = is_OFF(5); set_OFF(5) = false; % allow delete roi

    elseif any(popup_pairtype == keys_inspect_type_idx.nonmatched_ROIs) % NON-MATCHED ROIS 
        
        is_caller = popup_pairtype == keys_inspect_type_idx.nonmatched_ROIs;
        set_ON(2) = is_OFF(2) && roi_is_selected(~is_caller); set_OFF(2) = ~is_OFF(2) && ~roi_is_selected(~is_caller); % allow create pair
        set_ON(3) = is_OFF(3); set_OFF(3) = false; % allow out_of_fov roi
        set_ON(5) = is_OFF(5); set_OFF(5) = false; % allow delete roi
        
        
    elseif any(popup_pairtype == keys_inspect_type_idx.outofFOV_ROIs) % EDGES OR OUT OF FOVS ROIS
        set_ON(4) = is_OFF(4); set_OFF(4) = false; % allow add roi
        set_ON(5) = is_OFF(5); set_OFF(5) = false; % allow delete roi
        
    else % DELETED ROIS
        set_ON(4) = is_OFF(4); set_OFF(4) = false; % allow add roi
        set_ON(3) = is_OFF(3); set_OFF(3) = false; % allow out_of_fov roi
    end
    
end

% Set to OFF
if any(set_OFF)
    cellfun(@(h) set(h, 'Enable', 'off'), list_h(set_OFF))
end

% Set to ON
if any(set_ON)
    cellfun(@(h) set(h, 'Enable', 'on'), list_h(set_ON))
end


%__________________________________________________________________________
%
%   EDIT REG_RES - PUSHBUTTONS
%__________________________________________________________________________

% --- Executes on button press in pushbutton_match.
function pushbutton_match_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_match (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% NOTE: can match from: nonmatched_ROIs and missed_pairs

match_or_delete_pair(handles);


% --- Executes on button press in pushbutton_deletepair.
function pushbutton_deletepair_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deletepair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% NOTE: can delete pair from: matched_pairs, conflicting_pairs and
% incomplete_pairs

match_or_delete_pair(handles);



function match_or_delete_pair(handles)

% Get the type of inspection:
f = handles.keys.inspect_type_name{handles.gui_state.popup_pairtype};

% Get the rois to pair
pointer_rois = [get(handles.listbox_rois_s1, 'Value'), get(handles.listbox_rois_s2, 'Value')];
rois = [handles.keys.rois_pool{1}(pointer_rois(1)), handles.keys.rois_pool{2}(pointer_rois(2))];

% We need to sort the data according to the order in the current
% registration, not the order show on the GUI
rois = rois(handles.keys.sessions_reg_idx);
idx_session_pair = handles.keys.sessions_data_idx(handles.keys.sessions_reg_idx);

% Call special function to update REG_RES across multiple sessions
[handles.REG_RES] = gui_update_REG_RES_across_session({handles.REG_DATA.session}, handles.REG_RES, handles.keys.current_registration, idx_session_pair, rois, f);
% The pair might be moved to matched_pairs, conflicting_chain, or
% incomplete_chain. All affected pairs/rois will be moved accordingly.

% Re-Initialise rois tags (matched, unmatched, deleted, etc...)
handles.keys.rois_tag = init_ROIs_tags(handles.REG_RES(handles.keys.current_registration), handles.keys.sessions_reg_idx);

% Update the masks on tiffs
swap_multi_masks_rois_type([handles.uiaxes1, handles.uiaxes2], 1:handles.keys.n1n2(1), handles.keys.rois_tag{1}, handles.palette.tiff.linespec);
swap_multi_masks_rois_type([handles.uiaxes1, handles.uiaxes2], handles.keys.n1n2(1)+1:sum(handles.keys.n1n2), handles.keys.rois_tag{2}, handles.palette.tiff.linespec);

% Update pointer if inspect list of nonmatched_ROIs
if any( handles.gui_state.popup_pairtype == handles.keys.inspect_type_idx.nonmatched_ROIs ) 
    % Get caller
    is_caller = handles.gui_state.popup_pairtype == ...
        handles.keys.inspect_type_idx.(f);
    % Reset the non-ref to 1
    pointer_rois(~is_caller) = 1;
end

% Reset roi list and the selection 
[handles.keys.rois_pool, handles.keys.selected_rois] = get_rois_data_for_inspection(handles, [], handles.gui_state.popup_pairtype, 'pointer', pointer_rois);

% Update handles structure
guidata(handles.uifigure_regrois, handles);


% --- Executes on button press in pushbutton_addroi.
function pushbutton_addroi_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% NOTE: can add roi from deleted_ROIs and out_of_fov_ROIs

% Get the type of inspection:
f = handles.keys.inspect_type_name{handles.gui_state.popup_pairtype};

% Get caller 
is_caller = handles.gui_state.popup_pairtype == ...
    handles.keys.inspect_type_idx.(f);
% Identify where which session it is in the current registration
is_caller_ses_reg = is_caller(handles.keys.sessions_reg_idx);

% Get the roi to move
h_listROIs = {handles.listbox_rois_s1, handles.listbox_rois_s2};
pointer_roi_S0 = get(h_listROIs{is_caller}, 'Value');
roi_S0 = handles.keys.rois_pool{is_caller}(pointer_roi_S0);

% We need to sort the data according to the order in the current
% registration, not the order show on the GUI
rois = cell(1,2); 
rois{is_caller_ses_reg} = roi_S0;
idx_session_pair = handles.keys.sessions_data_idx(handles.keys.sessions_reg_idx);

% Call special function to update REG_RES across multiple sessions
[handles.REG_RES] = gui_update_REG_RES_across_session({handles.REG_DATA.session}, handles.REG_RES, handles.keys.current_registration, idx_session_pair, rois, f);
% The roi might be linked to matched_pairs, conflicting_chain, or
% incomplete_chain. All affected pairs/rois will be moved accordingly.

% Re-Initialise rois tags (matched, unmatched, deleted, etc...)
handles.keys.rois_tag = init_ROIs_tags(handles.REG_RES(handles.keys.current_registration), handles.keys.sessions_reg_idx);

% Update the masks on tiffs
swap_multi_masks_rois_type([handles.uiaxes1, handles.uiaxes2], 1:handles.keys.n1n2(1), handles.keys.rois_tag{1}, handles.palette.tiff.linespec);
swap_multi_masks_rois_type([handles.uiaxes1, handles.uiaxes2], handles.keys.n1n2(1)+1:sum(handles.keys.n1n2), handles.keys.rois_tag{2}, handles.palette.tiff.linespec);

% Update pointer 
pointer_rois = [1,1];
pointer_rois(is_caller) = pointer_roi_S0;

% Reset roi list and the selection 
[handles.keys.rois_pool, handles.keys.selected_rois] = get_rois_data_for_inspection(handles, [], handles.gui_state.popup_pairtype, 'pointer', pointer_rois);

% Update handles structure
guidata(handles.uifigure_regrois, handles);




% --- Executes on button press in pushbutton_deleteroi.
function pushbutton_deleteroi_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deleteroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% NOTE: we can delete a roi in 2 cases: if inspect nonmatched rois or if
% inspect out of fov rois

% New type to move the roi to
f_new = 'deleted_ROIs';


if handles.gui_state.popup_pairtype==handles.keys.inspect_type_idx.missed_pairs
    remove_2rois(handles, f_new)
else
    new_tag = -4;
    remove_single_roi(handles, f_new, new_tag)
end


% --- Executes on button press in pushbutton_outoffov.
function pushbutton_outoffov_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_outoffov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% NOTE: we can move to out of fov in 2 cases: if inspect nonmatched rois or
% if inspect deleted rois

% New type to move the roi to
f_new = 'outofFOV_ROIs';

if handles.gui_state.popup_pairtype==handles.keys.inspect_type_idx.missed_pairs
    remove_2rois(handles, f_new)
else
    new_tag = -3;
    remove_single_roi(handles, f_new, new_tag)
end



function remove_2rois(handles, f_new)

% Get the type of inspection:
f = handles.keys.inspect_type_name{handles.gui_state.popup_pairtype};

% Get the rois to pair
pointer_rois = [get(handles.listbox_rois_s1, 'Value'), get(handles.listbox_rois_s2, 'Value')];
rois = [handles.keys.rois_pool{1}(pointer_rois(1)), handles.keys.rois_pool{2}(pointer_rois(2))];

% We need to sort the data according to the order in the current
% registration, not the order show on the GUI
rois = rois(handles.keys.sessions_reg_idx);
idx_session_pair = handles.keys.sessions_data_idx(handles.keys.sessions_reg_idx);

% Call special function to update REG_RES across multiple sessions
[handles.REG_RES] = gui_update_REG_RES_across_session({handles.REG_DATA.session}, handles.REG_RES, handles.keys.current_registration, idx_session_pair, rois, f, f_new);
% The pair was moved to deleted rois/ out of fov. All affected pairs/rois
% will be moved accordingly.

% Re-Initialise rois tags (matched, unmatched, deleted, etc...)
handles.keys.rois_tag = init_ROIs_tags(handles.REG_RES(handles.keys.current_registration), handles.keys.sessions_reg_idx);

% Update the masks on tiffs
swap_multi_masks_rois_type([handles.uiaxes1, handles.uiaxes2], 1:handles.keys.n1n2(1), handles.keys.rois_tag{1}, handles.palette.tiff.linespec);
swap_multi_masks_rois_type([handles.uiaxes1, handles.uiaxes2], handles.keys.n1n2(1)+1:sum(handles.keys.n1n2), handles.keys.rois_tag{2}, handles.palette.tiff.linespec);

% Update pointer if inspect list of nonmatched_ROIs
if any( handles.gui_state.popup_pairtype == handles.keys.inspect_type_idx.nonmatched_ROIs ) 
    % Get caller
    is_caller = handles.gui_state.popup_pairtype == ...
        handles.keys.inspect_type_idx.(f);
    % Reset the non-ref to 1
    pointer_rois(~is_caller) = 1;
end

% Reset roi list and the selection 
[handles.keys.rois_pool, handles.keys.selected_rois] = get_rois_data_for_inspection(handles, [], handles.gui_state.popup_pairtype, 'pointer', pointer_rois);

% Update handles structure
guidata(handles.uifigure_regrois, handles);



function remove_single_roi(handles, f_new, new_tag)
% f_new : new type to move the roi to

% Get the type of inspection:
f = handles.keys.inspect_type_name{handles.gui_state.popup_pairtype};

% Get caller 
is_caller = handles.gui_state.popup_pairtype == ...
    handles.keys.inspect_type_idx.(f);
% Identify where which session it is in the current registration
is_caller_ses_reg = is_caller(handles.keys.sessions_reg_idx);

% Get the roi to move
h_listROIs = {handles.listbox_rois_s1, handles.listbox_rois_s2};
pointer_roi_S0 = get(h_listROIs{is_caller}, 'Value');
roi_S0 = handles.keys.rois_pool{is_caller}(pointer_roi_S0);

% Remove from original list 
is_roi_2change = handles.REG_RES(handles.keys.current_registration).(f){is_caller_ses_reg} == roi_S0;
handles.REG_RES(handles.keys.current_registration).(f){is_caller_ses_reg} = handles.REG_RES(handles.keys.current_registration).(f){is_caller_ses_reg}(~is_roi_2change);
% Add it to corresponding out_of_fov list
handles.REG_RES(handles.keys.current_registration).(f_new){is_caller_ses_reg} = ...
    [handles.REG_RES(handles.keys.current_registration).(f_new){is_caller_ses_reg}, roi_S0];

% Update the pointer of the caller
pointer_rois = [1,1];
pointer_rois(is_caller) = pointer_roi_S0;

% Change rois tag
handles.keys.rois_tag{is_caller}(roi_S0) = new_tag; 

% Update the masks on tiffs
ls_idx = roistag2linespec(new_tag);
if is_caller(1) % If inspect list of first session
    n_update = 0;
else
    n_update = handles.keys.n1n2(1);
end
swap_single_mask_roi_type([handles.uiaxes1, handles.uiaxes2], roi_S0+n_update, handles.palette.tiff.linespec(ls_idx))

% Reset roi list and the selection 
[handles.keys.rois_pool, handles.keys.selected_rois] = get_rois_data_for_inspection(handles, [], handles.gui_state.popup_pairtype, 'pointer', pointer_rois);

% Update handles structure
guidata(handles.uifigure_regrois, handles);



%__________________________________________________________________________
%
%   CALLBACKS VIEW PANEL
%__________________________________________________________________________
        


% --- Executes on button press in pushbutton_resetview.
function pushbutton_resetview_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_resetview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Simply reset the xy limits of each axis

% Only need to reset one axis if axes are linked, Otherwise reset the
% original limits of each axis
if handles.checkbox_linkaxes.Value == 0
    set(handles.uiaxes1, 'XLim', handles.orig_siz.plotlimits{1}.Xlim, 'YLim', handles.orig_siz.plotlimits{1}.Ylim)
    set(handles.uiaxes2, 'XLim', handles.orig_siz.plotlimits{2}.Xlim, 'YLim', handles.orig_siz.plotlimits{2}.Ylim)
else
    set(handles.uiaxes1, 'XLim', handles.orig_siz.plotlimits_max.Xlim, 'YLim', handles.orig_siz.plotlimits_max.Ylim)
end

% --- Executes on button press in togglebutton_hidetiff1.
function togglebutton_hidetiff1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_hidetiff1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_hidetiff1

if get(hObject,'Value')==0
    set(handles.uiaxes1.Children(end),'Visible', 'on');
    set(handles.uiaxes1.Children(end-2),'Visible', 'off');
else
    set(handles.uiaxes1.Children(end),'Visible', 'off');
    set(handles.uiaxes1.Children(end-2),'Visible', 'on');
end

% --- Executes on button press in togglebutton_hidetiff2.
function togglebutton_hidetiff2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_hidetiff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_hidetiff2

if get(hObject,'Value')==0
    set(handles.uiaxes2.Children(end),'Visible', 'on');
    set(handles.uiaxes2.Children(end-1),'Visible', 'off');
else
    set(handles.uiaxes2.Children(end),'Visible', 'off');
    set(handles.uiaxes2.Children(end-1),'Visible', 'on');
end

% --- Executes on button press in checkbox_linkaxes.
function checkbox_linkaxes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_linkaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_linkaxes

axes_are_linked = get(hObject,'Value')==1;

[img, imgoutline, rois] = get_session_data(handles, axes_are_linked);

% Change the tiff images
handles.uiaxes1.Children(end).CData = img{1};
handles.uiaxes2.Children(end).CData = img{2};

% Change the masks on tiff images
change_ROIs_masks(handles.uiaxes1.Children(1:end-3),[rois{1}{:}])
change_ROIs_masks(handles.uiaxes2.Children(1:end-3),[rois{2}{:}])

% Change the outlines image for the other session
change_image_outline(handles.uiaxes1.Children(end-1),imgoutline{1}{2})
change_image_outline(handles.uiaxes1.Children(end-2),imgoutline{1}{1})
change_image_outline(handles.uiaxes2.Children(end-1),imgoutline{2}{2})
change_image_outline(handles.uiaxes2.Children(end-2),imgoutline{2}{1})


if axes_are_linked
    
    % Link the axes
    linkaxes([handles.uiaxes1, handles.uiaxes2])
    
    % Re-enable 'mask aligned to'
    set([handles.radiobutton_s1, handles.radiobutton_s2], 'Enable', 'on')
    
else
    % de-link all axes
    linkaxes([handles.uiaxes1, handles.uiaxes2], 'off')
    % Disable 'mask aligned to'
    set([handles.radiobutton_s1, handles.radiobutton_s2], 'Enable', 'off')
end


function swap_ref_session(handles)
% This is only possible if axes are linked

% Swap ref
handles.keys.is_ref = ~handles.keys.is_ref;

% Update handles structure
guidata(handles.uifigure_regrois, handles);

% We can trigger 'radiobutton_linkaxes_Callback' to change the tiff data
% automatically.
checkbox_linkaxes_Callback(handles.checkbox_linkaxes, [], handles)


% --- Executes on button press in radiobutton_s1.
function radiobutton_s1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_s1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_s1

% Only do something if it was false
if ~handles.keys.is_ref(1)
    swap_ref_session(handles)
end

% --- Executes on button press in radiobutton_s2.
function radiobutton_s2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_s2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_s2

% Only do something if it was false
if ~handles.keys.is_ref(2)
    swap_ref_session(handles)
end




%__________________________________________________________________________
% 
%   CALLBACKS CHECKBOX SHOW/HIDE ROIS \
%        - CHOOSE TYPE TO DISPLAY: PAIRED, NON-MATCHED, OUT_OF_FOV, DELETED
%        - SHOW BOTH SESSIONS ON TIFF
%__________________________________________________________________________


% --- Executes on button press in checkbox_showboth.
function checkbox_showboth_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showboth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showboth

% Keep other session children
list_chd = [handles.uiaxes1.Children(handles.keys.n1n2(1)+1:end-3) ; handles.uiaxes2.Children(1:handles.keys.n1n2(1))];
list_chd_outline = [handles.uiaxes1.Children(end-1) ; handles.uiaxes2.Children(end-2)];
    
if hObject.Value==1
    on_or_off = 'on';
    
    % Check the rois types currently shown
    toogle_type = find([handles.checkbox_showmatched.Value==1, ...
        handles.checkbox_shownonmatched.Value==1, ...
        handles.checkbox_showoutoffov.Value==1, ...
        handles.checkbox_showdeleted.Value==1]);
   
    if isempty(toogle_type)
        list_chd_all = [];
    else
        % Find the rois with the corresponding tag
        is_roi_type_s1 = ismember(handles.keys.rois_tag{1}, handles.keys.tags_def{toogle_type(1)});
        is_roi_type_s2 = ismember(handles.keys.rois_tag{2}, handles.keys.tags_def{toogle_type(1)});
        for k=2:length(toogle_type)
            is_roi_type_s1 = is_roi_type_s1 | ismember(handles.keys.rois_tag{1}, handles.keys.tags_def{toogle_type(k)});
            is_roi_type_s2 = is_roi_type_s2 | ismember(handles.keys.rois_tag{2}, handles.keys.tags_def{toogle_type(k)});
        end
        
        list_is_type = [is_roi_type_s2' ; is_roi_type_s1'];

        list_chd_all = [list_chd(list_is_type);list_chd_outline];
    end
    
else
    on_or_off = 'off';
    list_chd_all = [list_chd;list_chd_outline];
end

if ~isempty(list_chd_all)
    arrayfun(@(h_chd) set(h_chd,'Visible', on_or_off), list_chd_all);
end


function toogle_rois_visibility(toogle_type, checkbox_value, chd_axes1, chd_axes2, tags_def, rois_tag, n1n2, checkbox_bothmasks_val)

if checkbox_value==1
    on_or_off = 'on';
else
    on_or_off = 'off';
end

% Find the rois with the tag
is_roi_type_s1 = ismember(rois_tag{1}', tags_def{toogle_type});
is_roi_type_s2 = ismember(rois_tag{2}', tags_def{toogle_type});

if checkbox_bothmasks_val==1
    % Keep all children
    list_chd = [chd_axes1(1:end-3); chd_axes2(1:end-3)];
    list_is_type = [is_roi_type_s1 ; is_roi_type_s2 ; is_roi_type_s1 ; is_roi_type_s2];  
else
    % Keep session only
    list_chd = [chd_axes1(1:n1n2(1)); chd_axes2(n1n2(1)+1:end-3)];
    list_is_type = [is_roi_type_s1 ; is_roi_type_s2];
end
list_chd = list_chd(list_is_type);

arrayfun(@(h_chd) set(h_chd,'Visible', on_or_off), list_chd);


% --- Executes on button press in checkbox_showmatched.
function checkbox_showmatched_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showmatched (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showmatched

toogle_rois_visibility(1, hObject.Value, handles.uiaxes1.Children, handles.uiaxes2.Children, ...
    handles.keys.tags_def, handles.keys.rois_tag, handles.keys.n1n2, handles.checkbox_showboth.Value)

% --- Executes on button press in checkbox_shownonmatched.
function checkbox_shownonmatched_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_shownonmatched (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_shownonmatched

toogle_rois_visibility(2, hObject.Value, handles.uiaxes1.Children, handles.uiaxes2.Children, ...
    handles.keys.tags_def, handles.keys.rois_tag, handles.keys.n1n2, handles.checkbox_showboth.Value)

% --- Executes on button press in checkbox_showoutoffov.
function checkbox_showoutoffov_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showoutoffov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showoutoffov

toogle_rois_visibility(3, hObject.Value, handles.uiaxes1.Children, handles.uiaxes2.Children, ...
    handles.keys.tags_def, handles.keys.rois_tag, handles.keys.n1n2, handles.checkbox_showboth.Value)

% --- Executes on button press in checkbox_showdeleted.
function checkbox_showdeleted_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showdeleted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showdeleted

toogle_rois_visibility(4, hObject.Value, handles.uiaxes1.Children, handles.uiaxes2.Children, ...
    handles.keys.tags_def, handles.keys.rois_tag, handles.keys.n1n2, handles.checkbox_showboth.Value)







% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%           VARIOUS GUI STUFF (CreateFcn...) NOT EDITED
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% --- Executes during object creation, after setting all properties.
function popupmenu_s1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_s1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_s2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_s2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_inspect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_inspect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function listbox_rois_s2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_rois_s2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function listbox_rois_s1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_rois_s1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





