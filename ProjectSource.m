classdef finalProject < matlab.apps.AppBase
 %just change the file name to the corresponding class name to fix the error
 % OR replace the class name with ProjectSource both should work
    % Properties that correspond to app components
    properties (Access = public)
        ImageprocessingLABbyBraveCakeUIFigure  matlab.ui.Figure
        PickanImageButton            matlab.ui.control.Button
        originalImage                matlab.ui.control.Image
        ProcessingtypeDropDownLabel  matlab.ui.control.Label
        filtersList                  matlab.ui.control.DropDown
        processButton                matlab.ui.control.Button
        resetSelection               matlab.ui.control.CheckBox
        ViewHistogramButton          matlab.ui.control.Button
        ViewFourierButton            matlab.ui.control.Button
        SaveButton                   matlab.ui.control.Button
    end

    
    properties (Access = public)
        imagePath % Description
    end
    
    properties (Access = private)
        result % Description
        originalSrc % Description
    end
    
    methods (Access = private)
        
        function results = rgb2grayScale(app,img)
            results = img;
            if(size(img,3)==3)
                results = rgb2gray(img);
            end
        end
        
        function results = getValidImg(app,path)
            [img,map] = imread(path);
            if(size(img,3)~=3)
                img = ind2rgb(img, map);
            end
            results = img;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            
        end

        % Button pushed function: PickanImageButton
        function PickanImageButtonPushed(app, event)
            [name,path] = uigetfile("*.jpg;*.png");
            if(name==0) %user selected nothing
                return
            end
            app.imagePath= [path name];
            app.originalImage.ImageSource = app.imagePath;
            app.result = getValidImg(app,app.imagePath);
            app.originalSrc=app.result; 
            % disp(app.originalImage.ImageSource); %display directory
            app.originalImage.Visible= 'on';
            app.processButton.Visible='on';
            app.resetSelection.Visible='on';
            app.ViewHistogramButton.Visible='on'; 
            app.ViewFourierButton.Visible='on';
            app.SaveButton.Visible='on';
        end

        % Button pushed function: processButton
        function processButtonPushed(app, event)
            if(app.resetSelection.Value)
                app.result = imread(app.originalImage.ImageSource); 
            end
            switch app.filtersList.Value
                case 'rgb2gray'
                    app.result= rgb2grayScale(app,app.result)
                case 'Scale up'
                    prompt = {'Height factor:','width factor:'};
                    hf_wf = inputdlg(prompt,"up sampling");
                    hf_wf= str2num(cell2mat(hf_wf));
                    if  (any(isnan(hf_wf)))
                        error('Invalid input: height and width factors must be numeric.')
                    end
                    app.result= Filters.up_sample(app.result,hf_wf(1),hf_wf(2));
                case 'Scale down'
                    prompt = {'Height factor:','width factor:'};
                    hf_wf= inputdlg(prompt,"down sampling");
                    hf_wf= str2num(cell2mat(hf_wf));
                    if  (any(isnan(hf_wf)))
                        error('Invalid input: height and width factors must be numeric.')
                    end
                    app.result= Filters.down_sample(app.result,hf_wf(1),hf_wf(2));
                case 'Negative Grayscale'
                    app.result = 255-rgb2grayScale(app,app.result);
                case 'Log intensity'
                    prompt = {'Enter the value of c:'};
                    c =  inputdlg(prompt,"Logarithmic intensity");
                    c= str2num(cell2mat(c))
                    if  (isnan(c)||isempty(c))
                        error('Invalid input c must be numeric.')
                    end
                    app.result = Filters.logIntensity(app.result,c(1));
                case 'Reverse log'
                    prompt = {'Enter the value of c:'};
                    c =  inputdlg(prompt,"Reverse logarithmic intensity");
                    if(isempty(c)||isempty(c))
                        return
                    end
                    c= str2num(cell2mat(c))
                    if  (isnan(c))
                        error('Invalid input c must be numeric.')
                    end
                    app.result = Filters.reverseLog(app.result,c(1)*10);
                case 'Power law'
                    prompt = {'Gamma:','C:'};
                    g_c= inputdlg(prompt,"Power law");
                    g_c= str2double(g_c);
                    if  (any(isnan(g_c))||isempty(g_c))
                        error('Invalid input it must be numeric.')
                    end
                    app.result= Filters.powerLawInesity(app.result,g_c(2),g_c(1));
                case 'Contrast Stretching'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'Minimum value:','Maximum value:'};
                    mn_mx= inputdlg(prompt,"Gray level Contrast Stetching");
                    mn_mx= str2double(mn_mx);
                    if  (any(isnan(mn_mx))||isempty(mn_mx))
                        error('Invalid input it must be numeric.')
                    end
                    app.result = Filters.conrastStretching(app.result,mn_mx(1),mn_mx(2))
                case 'Histogram Equalization'
                    app.result= Filters.histogramEqualization(rgb2grayScale(app,app.result));
                case 'Gray Bit Slicing'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'bit number:'};
                    bit= inputdlg(prompt,"bit slicing");
                    bit= str2num(bit{1});
                    if  (any(isnan(bit))||isempty(bit))
                        error('Invalid input it must be numeric.')
                    end
                    app.result= Filters.bitSlice(app.result,bit);
                case 'Graylevel Thresholding'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'Threshold pixel :'};
                    value= inputdlg(prompt,"thresholding");
                    value= str2num(value{1});
                    app.result = Filters.applyThresholding(app.result,value); %anything less than value is black
                case 'Gray Slicing'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'start pixel :','end pixel :'};
                    start_end= inputdlg(prompt,"graylevel slicing");
                    start_end= str2double(start_end);
                    [app.result,~] = Filters.graySlice(app.result,start_end(1),start_end(2));
                case 'Subtraction'
                    [name,path] = uigetfile("*.jpg;*.png");
                    if(name==0) %user selected nothing
                      return
                    end
                    subtrahendImg = getValidImg(app,[path name]);
                    app.result = app.result- subtrahendImg;
                case 'Addition'
                    [name,path] = uigetfile("*.jpg;*.png");
                    if(name==0) %user selected nothing
                      return
                    end
                    addedImg = getValidImg(app,[path name]);
                    app.result = app.result+ addedImg;
                case 'AND operation'
                    [name,path] = uigetfile("*.jpg;*.png");
                    if(name==0) %user selected nothing
                      return
                    end
                    addedImg = getValidImg(app,[path name]);
                    app.result = bitand(app.result, addedImg);
                case 'OR operation'
                    [name,path] = uigetfile("*.jpg;*.png");
                    if(name==0) %user selected nothing
                      return
                    end
                    addedImg = getValidImg(app,[path name]);
                    app.result = bitor(app.result,addedImg);
                case 'Average Filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'Neighbours number :'};
                    n= inputdlg(prompt,"Average Filter");
                    n= str2num(n{1});
                    app.result = Filters.averageFilter(app.result,n);
                case 'Median Filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'Neighbours number :'};
                    n= inputdlg(prompt,"Median Filter");
                    n= str2num(n{1});
                    app.result = Filters.medianFilter(app.result,n);
                case 'Max Filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'Neighbours number :'};
                    n= inputdlg(prompt,"Max Filter");
                    n= str2num(n{1});   
                    app.result = Filters.maxFilter(app.result,n);
                case 'Min Filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'Neighbours number :'};
                    n= inputdlg(prompt,"Min Filter");
                    n= str2num(n{1});   
                    app.result = Filters.minFilter(app.result,n);                     
                case 'Laplacian'
                    app.result = rgb2grayScale(app,app.result);
                    app.result = Filters.laplacian(app.result,false);
                case 'Composite Laplacian'
                    app.result = rgb2grayScale(app,app.result);
                    app.result = Filters.laplacian(app.result,true);
                case 'Robert'
                    app.result = rgb2grayScale(app,app.result);
                    app.result = Filters.robert(app.result);
                case 'Sobel'
                    app.result = rgb2grayScale(app,app.result);
                    app.result = Filters.sobel(app.result);
                case 'Low pass filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'D0 (cutoff):'};
                    d0= inputdlg(prompt,"ILPF");
                    d0= str2num(d0{1});   
                    app.result = Filters.ILPF(app.result,d0);
                case 'High pass filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'D0 (cutoff):'};
                    d0= inputdlg(prompt,"ILPF");
                    d0= str2num(d0{1});   
                    app.result = Filters.IHPF(app.result,d0);
                case 'Bitterwise low pass filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'D0 (cutoff):'};
                    d0= inputdlg(prompt,"ILPF");
                    d0= str2num(d0{1});   
                    app.result = Filters.BLPF(app.result,d0);
                case 'Bitterwise high pass filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'D0 (cutoff):'};
                    d0= inputdlg(prompt,"ILPF");
                    d0= str2num(d0{1});   
                    app.result = Filters.BHPF(app.result,d0);
                case 'Gaussian low pass filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'D0 (cutoff):'};
                    d0= inputdlg(prompt,"ILPF");
                    d0= str2num(d0{1});   
                    app.result = Filters.GLPF(app.result,d0);
                case 'Gaussian high pass filter'
                    app.result= rgb2grayScale(app,app.result);
                    prompt = {'D0 (cutoff):'};
                    d0= inputdlg(prompt,"ILPF");
                    d0= str2num(d0{1});   
                    app.result = Filters.GHPF(app.result,d0);
                    
                    
            end
            
            imshow(app.result);
        end

        % Button pushed function: ViewHistogramButton
        function ViewHistogramButtonPushed(app, event)
            Filters.compareHistogram(rgb2grayScale(app,app.originalSrc),rgb2grayScale(app,app.result))
        end

        % Button pushed function: ViewFourierButton
        function ViewFourierButtonPushed(app, event)
            Filters.compareFourier(rgb2grayScale(app,app.originalSrc),rgb2grayScale(app,app.result))
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            imwrite(app.result,'result.jpg');
            msgbox("Image saved successfully as result.jpg","saved successfully");
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ImageprocessingLABbyBraveCakeUIFigure and hide until all components are created
            app.ImageprocessingLABbyBraveCakeUIFigure = uifigure('Visible', 'off');
            app.ImageprocessingLABbyBraveCakeUIFigure.Position = [100 100 639 504];
            app.ImageprocessingLABbyBraveCakeUIFigure.Name = 'Image processing LAB ~ by BraveCake';

            % Create PickanImageButton
            app.PickanImageButton = uibutton(app.ImageprocessingLABbyBraveCakeUIFigure, 'push');
            app.PickanImageButton.ButtonPushedFcn = createCallbackFcn(app, @PickanImageButtonPushed, true);
            app.PickanImageButton.Position = [250 458 100 22];
            app.PickanImageButton.Text = 'Pick an Image';

            % Create originalImage
            app.originalImage = uiimage(app.ImageprocessingLABbyBraveCakeUIFigure);
            app.originalImage.Visible = 'off';
            app.originalImage.Position = [92 43 476 282];

            % Create ProcessingtypeDropDownLabel
            app.ProcessingtypeDropDownLabel = uilabel(app.ImageprocessingLABbyBraveCakeUIFigure);
            app.ProcessingtypeDropDownLabel.HorizontalAlignment = 'right';
            app.ProcessingtypeDropDownLabel.Position = [184 408 91 22];
            app.ProcessingtypeDropDownLabel.Text = 'Processing type';

            % Create filtersList
            app.filtersList = uidropdown(app.ImageprocessingLABbyBraveCakeUIFigure);
            app.filtersList.Items = {'rgb2gray', 'Scale down', 'Scale up', 'Negative Grayscale', 'Log intensity', 'Reverse log', 'Buffer', 'Power law', 'Contrast Stretching', 'Gray Bit Slicing', 'Histogram Equalization', 'Graylevel Thresholding', 'Gray Slicing', 'Subtraction', 'Addition', 'AND operation', 'OR operation', 'Average Filter', 'Median Filter', 'Max Filter', 'Min Filter', 'Laplacian', 'Composite Laplacian', 'Robert', 'Sobel', 'Low pass filter', 'High pass filter', 'Bitterwise low pass filter', 'Bitterwise high pass filter', 'Gaussian low pass filter', 'Gaussian high pass filter'};
            app.filtersList.Position = [290 408 100 22];
            app.filtersList.Value = 'rgb2gray';

            % Create processButton
            app.processButton = uibutton(app.ImageprocessingLABbyBraveCakeUIFigure, 'push');
            app.processButton.ButtonPushedFcn = createCallbackFcn(app, @processButtonPushed, true);
            app.processButton.Visible = 'off';
            app.processButton.Position = [250 356 100 22];
            app.processButton.Text = 'Process';

            % Create resetSelection
            app.resetSelection = uicheckbox(app.ImageprocessingLABbyBraveCakeUIFigure);
            app.resetSelection.Visible = 'off';
            app.resetSelection.Text = 'Reset all previous processing';
            app.resetSelection.Position = [417 408 179 22];

            % Create ViewHistogramButton
            app.ViewHistogramButton = uibutton(app.ImageprocessingLABbyBraveCakeUIFigure, 'push');
            app.ViewHistogramButton.ButtonPushedFcn = createCallbackFcn(app, @ViewHistogramButtonPushed, true);
            app.ViewHistogramButton.Visible = 'off';
            app.ViewHistogramButton.Position = [457 377 100 22];
            app.ViewHistogramButton.Text = 'View Histogram';

            % Create ViewFourierButton
            app.ViewFourierButton = uibutton(app.ImageprocessingLABbyBraveCakeUIFigure, 'push');
            app.ViewFourierButton.ButtonPushedFcn = createCallbackFcn(app, @ViewFourierButtonPushed, true);
            app.ViewFourierButton.Visible = 'off';
            app.ViewFourierButton.Position = [457 335 100 22];
            app.ViewFourierButton.Text = 'View Fourier';

            % Create SaveButton
            app.SaveButton = uibutton(app.ImageprocessingLABbyBraveCakeUIFigure, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.Visible = 'off';
            app.SaveButton.Position = [250 12 100 22];
            app.SaveButton.Text = 'Save';

            % Show the figure after all components are created
            app.ImageprocessingLABbyBraveCakeUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = finalProject

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.ImageprocessingLABbyBraveCakeUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ImageprocessingLABbyBraveCakeUIFigure)
        end
    end
end