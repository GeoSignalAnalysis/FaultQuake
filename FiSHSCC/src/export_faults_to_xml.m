function export_faults_to_xml(DATA)

    % Check if './Sources' folder exists, if not, create it
%     if ~exist('../Sources', 'dir')
%         mkdir('./Sources');
%     end

    % Loop over the faults
    fault_names = fieldnames(DATA);
    for i = 1:length(fault_names)
        fault = DATA.(fault_names{i});
        
        % Determine XML filename from fault name and save in './Sources' folder
        xml_filename = ['./output_files/Fault_OQ/' fault_names{i} '.xml'];
        currentDirectory = pwd;
        disp(['Current directory: ' currentDirectory]);
        
        % Open the XML file for writing
        fid = fopen(xml_filename, 'w');


    % Loop over the faults
 

        % Write the header information
        fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
        fprintf(fid, '<nrml xmlns="http://openquake.org/xmlns/nrml/0.4" xmlns:gml="http://www.opengis.net/gml">\n');
        fprintf(fid, '    <sourceModel name="%s">\n', fault_names{i});

        % Write the fault details to the XML
        fprintf(fid, '        <simpleFaultSource id="%d" name="Simple Fault Source" tectonicRegion="Active Shallow Crust">\n', fault.id);
        fprintf(fid, '            <simpleFaultGeometry>\n');
        fprintf(fid, '                <gml:LineString>\n');
        fprintf(fid, '                    <gml:posList>\n');
        fprintf(fid, '                        <!-- Fill in coordinates here -->\n');
        fprintf(fid, '                    </gml:posList>\n');
        fprintf(fid, '                </gml:LineString>\n');
        fprintf(fid, '                <dip>%d</dip>\n', fault.Dip);
        fprintf(fid, '                <upperSeismoDepth>%e</upperSeismoDepth>\n', fault.Seismogenic_Thickness - fault.Telap);
        fprintf(fid, '                <lowerSeismoDepth>%e</lowerSeismoDepth>\n', fault.Seismogenic_Thickness);
        fprintf(fid, '            </simpleFaultGeometry>\n');
        fprintf(fid, '            <magScaleRel>%s</magScaleRel>\n', fault.ScR);
        fprintf(fid, '            <ruptAspectRatio>%e</ruptAspectRatio>\n', 2.0000000E+00);  % This is a placeholder, adjust if necessary
%         fprintf(fid, '            <incrementalMFD minMag="%f" binWidth="%f">\n', fault.Mmin, fault.bin);
        fprintf(fid, '                <occurRates>%s</occurRates>\n', num2str(fault.rates, '%e '));
        fprintf(fid, '            </incrementalMFD>\n');
        fprintf(fid, '            <rake>%e</rake>\n', 9.0000000E+01);  % This is a placeholder, adjust if necessary
        fprintf(fid, '        </simpleFaultSource>\n');

        % Finish the XML
        fprintf(fid, '    </sourceModel>\n');
        fprintf(fid, '</nrml>\n');

        % Close the XML file
        fclose(fid);
    end
end

