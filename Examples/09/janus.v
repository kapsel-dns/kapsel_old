APPS.MultiWindowApp MultiWindowApp<NEscalingFactor=0.9561128616,NExOffset=-9.,NEyOffset=13.08687115> {
   UI {
      shell {
         x = 2220;
         y = 88;
         cancel = 1;
      };
      Modules {
         IUI {
            optionList {
               cmdList => {
                  <-.<-.<-.<-.Read_Field.read_field_ui.panel.option,
                  <-.<-.<-.<-.bounds.UIpanel.option,
                  <-.<-.<-.<-.Basic_Axis3D.panel.option,
                  <-.<-.<-.<-.isosurface.UIpanel.option,
                  <-.<-.<-.<-.downsize.panel.option,
                  <-.<-.<-.<-.combine_vect.CombVectUI.UIpanel.option,
                  <-.<-.<-.<-.glyph.UIpanel.option,
                  <-.<-.<-.<-.image_capture.UImod_panel.option,
                  <-.<-.<-.<-.isosurface#1.UIpanel.option,
                  <-.<-.<-.<-.extract_component.ExtrCompUI.panel.option,
                  <-.<-.<-.<-.data_math.UIpanel.option,
                  <-.<-.<-.<-.isosurface#2.UIpanel.option
               };
               selectedItem = 0;
            };
            mod_panel {
               x = 0;
               y = 0;
            };
         };
      };
      exit_dialog {
         ExitDialog {
            x = 1163;
            y = 625;
            width = 234;
            height = 112;
            cancel = 1;
         };
      };
   };
   GDM.Uviewer3D Uviewer3D<NEx=468.,NEy=738.> {
      Scene {
         Top {
            child_objs => {
               <-.<-.<-.Basic_Axis3D.out_obj,<-.<-.<-.bounds.out_obj,
               <-.<-.<-.glyph.out_obj,<-.<-.<-.isosurface.out_obj};
            AltObject<instanced=0>;
            Xform {
               ocenter = {27.2448,27.2448,
28.};
               dcenter = {-126.759,-89.7789,
-117.176};
               mat = {
                  0.0812481,-0.00302243,0.0304711,0.,0.0306207,0.00820233,-0.0808332,
0.,-6.51248e-05,0.0863857,0.0087419,0.,0.,0.,0.,1.
               };
               xlate = {126.757,89.7012,
117.36};
               center = {27.2448,27.2448,
28.};
            };
            Props {
               col = {0.,0.,0.};
            };
         };
         View {
            View {
               trigger = 1;
               back_col = {1.,1.,1.};
               refresh = 1;
            };
            ViewUI {
               FullViewWindow<instanced=0>;
               ViewPanel {
                  UI {
                     panel {
                        defaultHeight = 874;
                        defaultX = 705;
                        defaultY = 429;
                        defaultWidth = 977;
                     };
                  };
               };
            };
         };
      };
      Scene_Editor<NEdisplayMode="open"> {
         Track_Editor {
            GDxform_edit {
               abs_x_cent = 27.24475098;
               abs_y_cent = 27.24475098;
               abs_z_cent = 28.;
               abs_x_trans = 1.288341522;
               abs_y_trans = -16.88174629;
               abs_z_trans = 60.64237213;
               abs_x_rot = 26.89226723;
               abs_y_rot = -75.18087006;
               abs_z_rot = 117.179245;
               abs_scale = 0.1041243523;
               scale => ;
            };
            IUI {
               optionList {
                  selectedItem = 0;
               };
               Interactor_Behavior {
                  IUI {
                     panel {
                        visible = 1;
                     };
                     RotateMode {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           ,};
                     };
                     XYTranslateMode {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           ,,};
                     };
                  };
               };
            };
         };
         Camera_Editor {
            GDcamera_edit {
               front = -88.;
            };
            IUI {
               optionList {
                  selectedItem = 0;
               };
               General {
                  IUI {
                     panel {
                        visible = 1;
                     };
                     CameraType {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           };
                     };
                     Extents {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           };
                     };
                     Mapping {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           };
                     };
                     Normalize {
                        OPcmdList = {,,,
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           }};
                     };
                     DepthSort {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           ,};
                     };
                  };
               };
            };
         };
         Print_Editor {
            IUI<instanced=0> {
               Format {
                  OPcmdList = {
                     ,,
                     {
                        set=1,,,,,,,,,,,,,,,,,,,,
                     },
                     ,,
                  };
               };
               Orientation {
                  OPcmdList = {
                     {
                        set=1,,,,,,,,,,,,,,,,,,,,
                     },
                     };
               };
               Background {
                  OPcmdList = {
                     {
                        set=1,,,,,,,,,,,,,,,,,,,,
                     },
                     ,};
               };
               SizeMenu {
                  OPcmdList = {
                     ,
                     {
                        set=1,,,,,,,,,,,,,,,,,,,,
                     },
                     ,,,
                  };
               };
            };
         };
         View_Editor {
            IUI {
               optionList {
                  selectedItem = 0;
               };
               General {
                  IUI {
                     panel {
                        visible = 1;
                     };
                     Renderer {
                        OPcmdList = {
                           ,
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           ,,
                        };
                     };
                     Color {
                        Imm {
                           set = 1;
                        };
                        ColorEcho {
                           ColorView {
                              trigger = 1;
                           };
                        };
                        rgb_or_hsv {
                           s = 0.;
                           v = 1.;
                        };
                        x = 0;
                     };
                  };
               };
               Options {
                  IUI {
                     DoubleBuffer {
                        OPcmdList = {,,
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           }};
                     };
                     Aspect {
                        OPcmdList = {,
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           }};
                     };
                     BlendMode {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           };
                     };
                  };
               };
               Output {
                  IUI {
                     OutFBType {
                        OPcmdList = {,
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           ,};
                     };
                     OutZBType {
                        OPcmdList = {
                           ,,
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           ,
                        };
                     };
                  };
               };
            };
         };
         Light_Editor {
            IUI {
               optionList {
                  selectedItem = 0;
               };
               General {
                  IUI {
                     panel {
                        visible = 1;
                     };
                     AllFrame {
                        y = 0;
                     };
                     LightType {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           ,,};
                     };
                     VUIColorEditor {
                        ColorEcho {
                           ColorView {
                              trigger = 2;
                           };
                        };
                        x = 0;
                     };
                  };
               };
            };
         };
         Object_Editor {
            GDprops_edit {
               inherit = 0;
            };
            IUI {
               optionList {
                  selectedItem = 0;
               };
               General {
                  IUI {
                     panel {
                        visible = 1;
                     };
                     Pickable {
                        OPcmdList = {
                           ,
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           ,,
                        };
                     };
                     AltSpace {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           };
                     };
                     TransformMode {
                        OPcmdList = {
                           {
                              set=1,,,,,,,,,,,,,,,,,,,,
                           },
                           ,,};
                     };
                  };
               };
               Properties {
                  IUI {
                     ObjectOptions = {
                        {
                           set=1,,,,,,,,,,,,,,,,,,,,
                        },
                        };
                     EditAltProps {
                        y = 0;
                        x = 0;
                     };
                     Type {
                        IUI {
                           optionList {
                              selectedItem = 0;
                           };
                           General {
                              IUI {
                                 panel {
                                    visible = 1;
                                 };
                                 ColorOptions = {
                                    {
                                       set=1,,,,,,,,,,,,,,,,,,,,
                                    },
                                    ,};
                                 ColorEditor {
                                    Imm {
                                       set = 1;
                                    };
                                    ColorEcho {
                                       ColorView {
                                          trigger = 2;
                                       };
                                    };
                                    rgb_or_hsv {
                                       v = 0.;
                                    };
                                    x = 0;
                                 };
                              };
                           };
                        };
                     };
                  };
               };
            };
         };
         Datamap_Editor {
            IUI {
               Options {
                  IUI {
                     optionList {
                        selectedItem = 0;
                     };
                     Edit_Color {
                        IUI {
                           panel {
                              visible = 1;
                           };
                           ColorRB {
                              OPcmdList = {
                                 {
                                    set=1,,,,,,,,,,,,,,,,,,,,
                                 },
                                 };
                           };
                        };
                     };
                     Edit_Alpha {
                        IUI {
                           AlphaRB {
                              OPcmdList = {,
                                 {
                                    set=1,,,,,,,,,,,,,,,,,,,,
                                 }};
                           };
                        };
                     };
                     Edit_Range_Data {
                        IUI {
                           EditMenu {
                              OPcmdList = {
                                 {
                                    set=1,,,,,,,,,,,,,,,,,,,,
                                 },
                                 };
                           };
                        };
                     };
                     InputOutput {
                        IUI {
                           FileOP {
                              OPcmdList = {
                                 {
                                    set=1,,,,,,,,,,,,,,,,,,,,
                                 },
                                 };
                           };
                        };
                     };
                  };
               };
               DmapEcho {
                  DmapView {
                     trigger = 2;
                  };
               };
               ModelOptions = {
                  {
                     set=1,,,,,,,,,,,,,,,,,,,,
                  },
                  };
            };
         };
         Graph_Editor<NEdisplayMode="open"> {
            IUI {
               Properties<NEdisplayMode="open">;
            };
         };
      };
   };
   MODS.Read_Field Read_Field<NEx=9.,NEy=36.> {
      read_field_ui {
         file_browser {
            x = 819;
            y = 800;
            width = 300;
            height = 386;
            ok = 1;
            dirMaskCache = "/home/john/Work/KAPSEL/self_propulsion/swimmer_visualization/N6_single_a+2.0/avs_ch/*";
            cancel = 1;
         };
         portable = 0;
         filename = "/home/john/Work/KAPSEL/self_propulsion/swimmer_visualization/N6_single_a+2.0/avs_ch/data.fld";
         panel {
            option {
               set = 1;
            };
         };
      };
      DVread_field<NEdisplayMode="open"> {
         Mesh_Unif+Node_Data Output_Field;
         stepno<NEdisplayMode="open">;
      };
      Read_Field_Param<NEdisplayMode="open"> {
         current_step<NEportLevels={0,3}> = 1;
         one_time = 0;
         continuous = 0;
         bounce = 0;
      };
      DataObject {
         AltObject<instanced=0>;
      };
      loop {
         done = 1;
      };
      do_loop {
         set_run {
            output = 0;
         };
      };
   };
   MODS.bounds bounds<NEx=45.,NEy=630.> {
      in_field => <-.Read_Field.field;
      obj {
         AltObject<instanced=0>;
      };
      BoundsParam {
         imin = 1;
         imax = 1;
         jmin = 1;
         jmax = 1;
         kmin = 1;
         kmax = 1;
      };
   };
   HLM.Basic_Axis3D Basic_Axis3D<NEx=432.,NEy=864.> {
      axis {
         xform {
            xlate = {-7.,-7.,-7.};
            mat = {
               10.07,0.,0.,0.,0.,10.07,0.,0.,0.,0.,10.07,0.,0.,0.,0.,1.
            };
         };
      };
      GroupObject {
         AltObject<instanced=0>;
      };
      text_glyph {
         obj {
            AltObject<instanced=0>;
         };
      };
      Axis_UI<instanced=0> {
         probe_edit {
            GDxform_editor {
               x_trans = -7.;
               y_trans = -7.;
               z_trans = -7.;
               abs_z_trans = -7.;
               abs_y_trans = -7.;
               abs_x_trans = -7.;
               scale = 10.07;
               abs_scale = 10.06999969;
            };
            XformEditorUI {
               trans_shell {
                  x = 853;
                  y = 480;
                  cancel = 1;
                  ok = 1;
               };
            };
         };
      };
      panel {
         option {
            set = 0;
         };
      };
   };
   MODS.isosurface isosurface<NEx=243.,NEy=396.> {
      in_field => <-.data_math.out_fld;
      IsoParam {
         iso_level => 0.56;
      };
      DVcell_data_labels {
         labels[];
      };
      DVnode_data_labels {
         labels[];
      };
      IsoUI {
         UIoptionBoxLabel {
            label_cmd {
               cmd[1];
            };
         };
         UIiso_level_typein {
            valEditor {
               y => 603;
               UI {
                  editor_shell {
                     x => 1523;
                  };
               };
            };
         };
      };
      obj {
         AltObject<instanced=0>;
         Datamap {
            DataRange = {
               {
                  DataMaxValue=2.,,DataMinValue=1.,,,,,,,,,
               }};
            DatamapValue = {
               {
                  value=1.,,v3=0.7581981421,v2=0.5804991722,
               },
               {
                  value=2.,v4=0.,v3=0.,,
               }};
         };
      };
   };
   MODS.downsize downsize<NEx=72.,NEy=315.> {
      in_field => <-.combine_vect.out_fld;
      DownsizeParam {
         factor0 = 2.;
         factor1 = 2.;
         factor2 = 3.;
      };
      obj {
         AltObject<instanced=0>;
      };
   };
   MODS.combine_vect combine_vect<NEx=72.,NEy=234.> {
      CombVectUI {
         DVnode_data_labels {
            labels[];
         };
         UIoptionBoxLabel {
            label_cmd {
               cmd[10];
            };
         };
      };
      obj {
         AltObject<instanced=0>;
         Datamap {
            DataRange = {
               {
                  DataMaxValue=255.,,DataMinValue=0.,,,,,,,,,
               }};
            DatamapValue = {
               {
                  value=0.,,,v2=0.1666666716,
               },
               {
                  value=255.,,,,
               }};
         };
      };
      in_field => <-.Read_Field.field;
   };
   MODS.glyph glyph<NEx=72.,NEy=396.> {
      in_field => <-.downsize.out_fld;
      in_glyph => <-.Arrow1.out_fld;
      GlyphParam {
         scale = 250.;
      };
      GlyphUI {
         scale_slider {
            min = 0.;
            max = 500.;
         };
         scale_slider_typein {
            valEditor {
               y => 591;
               UI {
                  editor_shell {
                     x => 359;
                  };
               };
            };
         };
      };
      obj {
         AltObject<instanced=0>;
         Datamap {
            DataRange[2] = {
               {
                  DataMaxValue=0.002259901259,,,,,,,,,,,,,,DataMinValue=2.260397923e-05,,,,,,,,,
               },
               {
                  DataMaxValue=0.02239557914,,DataMinValue=><-.DataRange[0].DataMaxValue,controlPoints=>
                  {DatamapValue[2],
                     DatamapValue[3]},,,,,,,,
               }};
            DatamapValue[4] = {
               {
                  value=2.260397923e-05,,v3=0.04079173505,v2=0.7097603083,v1=1.
               },
               {
                  value=0.002259901259,,,,v1=0.
               },
               {
                  value=0.002259901259,v4=1.,v3=1.,v2=0.,v1=0.
               },
               {
                  value=0.02239557914,v4=1.,v3=1.,v2=0.,v1=0.
               }};
         };
      };
   };
   GEOMS.Arrow1 Arrow1<NEx=72.,NEy=513.>;
   ANIM_MODS.image_capture image_capture<NEx=666.,NEy=792.> {
      imcapUI {
         MovieControls {
            ChooseDefaults {
               set = 0;
            };
            BR_QTY_Menu {
               selectedItem = 1;
            };
            Bitrate {
               value = 8.;
            };
            Quality {
               value = 23.;
            };
            Buffer {
               value = 69.;
            };
            FileSelector {
               fileFB {
                  width = 313;
                  height = 386;
                  y = 67;
                  ok = 1;
                  dirMaskCache = "/home/john/Work/KAPSEL/rigid_body/movie/*";
               };
            };
            ChooseBitrate {
               set = 1;
            };
         };
      };
      imcapCompute {
         ImageCap {
            LGDView => <-.<-.<-.Uviewer3D.Scene_Selector.curr_view;
         };
      };
      imcapParam {
         capture_mode = 0;
         tempdir = "/home/john/temp/";
         filename = "/home/john/temp/push_a0.2.mpg";
      };
   };
   MODS.isosurface isosurface#1<NEx=423.,NEy=234.> {
      in_field => <-.Read_Field.field;
      IsoParam {
         iso_component = 3;
         iso_level => -0.49;
         map_component = {3};
      };
      DVcell_data_labels {
         labels[];
      };
      DVnode_data_labels {
         labels[];
      };
      IsoUI {
         UIoptionBoxLabel {
            label_cmd {
               cmd[10];
            };
         };
         UIiso_level {
            min => -1.;
            decimalPoints = 2;
         };
         UIiso_level_typein {
            valEditor {
               y => 168;
            };
         };
      };
      obj {
         AltObject<instanced=0>;
         Datamap {
            DataRange[2] = {
               {
                  DataMaxValue=-0.8026258349,,,,,,,,,,,,,,DataMinValue=-1.,,,,,,,,,
               },
               {
                  DataMaxValue=1.,,DataMinValue=><-.DataRange[0].DataMaxValue,controlPoints=>
                  {DatamapValue[2],
                     DatamapValue[3]},,,,,,,,
               }};
            DatamapValue[4] = {
               {
                  value=-1.,,,v2=0.6666666865,v1=1.
               },
               {
                  value=-0.8026258349,,v3=0.1881539375,v2=0.2151977867,v1=0.
               },
               {
                  value=-0.8026258349,v4=1.,v3=0.1881539375,v2=0.2151977867,v1=0.
               },
               {
                  value=1.,v4=1.,v3=1.,v2=0.,v1=0.
               }};
         };
      };
      Iso {
         DVnmap {
            out {
               nnode_data = 1;
            };
         };
      };
   };
   MODS.extract_component extract_component<NEx=243.,NEy=234.> {
      in_field => <-.Read_Field.field;
      ExtrCompParam {
         component = 3;
      };
      ExtrCompUI {
         DVnode_data_labels {
            labels[];
         };
         UIradioBoxLabel {
            label_cmd {
               cmd[];
            };
         };
      };
      obj {
         AltObject<instanced=0>;
      };
   };
   MODS.data_math data_math<NEx=243.,NEy=315.> {
      in_field1 => <-.extract_component.out_fld;
      obj {
         AltObject<instanced=0>;
      };
      expres = "abs(#1)";
   };
   CMAP_EDTR.ColorMapEditor Janus_inside<NEx=243.,NEy=459.> {
      InObj => <-.isosurface.out_obj;
      TransferFunction {
         PointEditor {
            real_scalar = 1.;
         };
         Scene {
            Top {
               AltObject<instanced=0>;
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      UIshell {
         x = 1774;
         y = 358;
      };
      DmapRamp {
         Scene {
            Top {
               AltObject<instanced=0>;
               Xform {
                  ocenter = {170.,16.,0.};
                  dcenter = {-2957.,-251.601,0.};
               };
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      EditFields {
         UIoptionBoxRefresh {
            selectedItems = {0};
         };
         UIoptionRefresh {
            set = 1;
         };
      };
      ColorWheel {
         Scene {
            Top {
               AltObject<instanced=0>;
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
         HSPosConverter {
            hue = 0.5804991709;
            saturation = 0.7581981421;
         };
         PickField {
            xform {
               xlate = {-0.663265,-0.367347,
0.};
            };
         };
      };
      ValueCone {
         Scene {
            Top {
               AltObject<instanced=0>;
               Xform {
                  ocenter = {0.15,0.51,0.};
                  dcenter = {25.2413,85.82,
0.};
               };
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      MakeDmap {
         dmap_out {
            dataMax = 1.;
            DataRange = {
               {
                  DataMinValue=1.,,,,,DataMaxValue=2.,,,,,,,,,,,,controlPoints=>
                  {DatamapValue[0],
                     DatamapValue[1]},
               }};
            DatamapValue = {
               {
                  v1=0.,v2=0.5804991722,v3=0.7581981421,v4=1.,value=1.
               },
               {
                  v1=1.,v2=0.,v3=0.,v4=0.,value=2.
               }};
         };
         regen = 0;
         delta = 0.003921568859;
      };
      Presets {
         preset_load {
            preset = {0.,0.,1.,1.,
1.,1.,1.,0.,0.,0.};
         };
      };
   };
   CMAP_EDTR.ColorMapEditor Janus_bottom<NEx=423.,NEy=315.> {
      InObj => <-.isosurface#1.out_obj;
      TransferFunction {
         Scene {
            Top {
               AltObject<instanced=0>;
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
         PointEditor {
            real_scalar = -1.;
         };
      };
      UIshell {
         x = 939;
         y = 859;
         iconic = 1;
      };
      DmapRamp {
         Scene {
            Top {
               AltObject<instanced=0>;
               Xform {
                  ocenter = {170.,16.,0.};
                  dcenter = {-2957.,-251.601,0.};
               };
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      EditFields {
         UIoptionBoxRefresh {
            selectedItems = {0};
         };
         UIoptionRefresh {
            set = 1;
         };
      };
      ColorWheel {
         Scene {
            Top {
               AltObject<instanced=0>;
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
         HSPosConverter {
            hue = 0.6666666667;
            saturation = 1.;
         };
         PickField {
            xform {
               xlate = {-0.5,-0.866025,0.};
            };
         };
      };
      ValueCone {
         Scene {
            Top {
               AltObject<instanced=0>;
               Xform {
                  ocenter = {0.15,0.51,0.};
                  dcenter = {25.2413,85.82,
0.};
               };
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      MakeDmap {
         dmap_out {
            dataMax = 1.;
            dataMin = -1.;
            DataRange[2] = {
               {
                  DataMinValue=-1.,,,,,DataMaxValue=-0.8026258349,,,,,,,,,,,,controlPoints=>
                  {DatamapValue[0],
                     DatamapValue[1]},
               },
               {
                  DataMinValue=><-.DataRange[0].DataMaxValue,,,,,DataMaxValue=1.,,,,,,,,,,,,controlPoints=>
                  {DatamapValue[2],
                     DatamapValue[3]},
               }};
            DatamapValue[4] = {
               {
                  v1=1.,v2=0.6666666865,v3=1.,v4=1.,value=-1.
               },
               {
                  v1=0.,v2=0.2151977867,v3=0.1881539375,v4=1.,value=-0.8026258349
               },
               {
                  v1=0.,v2=0.2151977867,v3=0.1881539375,v4=1.,value=-0.8026258349
               },
               {
                  v1=0.,v2=0.,v3=1.,v4=1.,value=1.
               }};
         };
         regen = 0;
         delta = 0.007843137719;
      };
      Presets {
         preset_load {
            preset = {
               0.,1.,0.,0.,1.,0.1,0.,1.,0.,0.,1.,0.,1.,0.,0.
            };
         };
      };
   };
   MODS.isosurface isosurface#2<NEx=585.,NEy=234.> {
      in_field => <-.Read_Field.field;
      IsoParam {
         iso_component = 3;
         iso_level => 0.54;
         map_component = {3};
      };
      DVcell_data_labels {
         labels[];
      };
      DVnode_data_labels {
         labels[];
      };
      IsoUI {
         UIoptionBoxLabel {
            label_cmd {
               cmd[10];
            };
         };
         UIiso_level {
            min => -1.;
            decimalPoints = 2;
         };
         UIiso_level_typein {
            valEditor {
               y => 168;
            };
         };
      };
      obj {
         AltObject<instanced=0>;
         Datamap {
            DataRange[2] = {
               {
                  DataMaxValue=0.,,,,,,,,,,,,,,DataMinValue=-1.,,,,,,,,,
               },
               {
                  DataMaxValue=1.,,DataMinValue=><-.DataRange[0].DataMaxValue,controlPoints=>
                  {DatamapValue[2],
                     DatamapValue[3]},,,,,,,,
               }};
            DatamapValue[4] = {
               {
                  value=-1.,,v3=0.8271973729,v2=0.6348015666,
               },
               {
                  value=0.,,,v2=0.6666666865,v1=0.
               },
               {
                  value=0.,v4=1.,v3=1.,v2=0.,v1=1.
               },
               {
                  value=0.,v4=1.,v3=0.3479741216,v2=0.07732812315,v1=0.9887584448
               }};
         };
      };
      Iso {
         DVnmap {
            out {
               nnode_data = 1;
            };
         };
      };
   };
   CMAP_EDTR.ColorMapEditor Janus_top<NEx=585.,NEy=315.> {
      InObj => <-.isosurface#2.out_obj;
      TransferFunction {
         Scene {
            Top {
               AltObject<instanced=0>;
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
         PointEditor {
            real_scalar = -1.;
         };
      };
      UIshell {
         x = 1199;
         y = 333;
         iconic = 1;
      };
      DmapRamp {
         Scene {
            Top {
               AltObject<instanced=0>;
               Xform {
                  ocenter = {170.,16.,0.};
                  dcenter = {-3122.,-267.13,0.};
               };
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      EditFields {
         UIoptionBoxRefresh {
            selectedItems = {0};
         };
         UIoptionRefresh {
            set = 1;
         };
      };
      ColorWheel {
         Scene {
            Top {
               AltObject<instanced=0>;
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
         HSPosConverter {
            hue = 0.6348015467;
            saturation = 0.8271973729;
         };
         PickField {
            xform {
               xlate = {-0.547809,-0.619807,
0.};
            };
         };
      };
      ValueCone {
         Scene {
            Top {
               AltObject<instanced=0>;
               Xform {
                  ocenter = {0.15,0.51,0.};
                  dcenter = {26.5619,90.31,
0.};
               };
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      MakeDmap {
         dmap_out {
            dataMax = 1.;
            dataMin = -1.;
            DataRange[3] = {
               {
                  DataMinValue=-1.,,,,,DataMaxValue=0.,,,,,,,,,,,,controlPoints=>
                  {DatamapValue[0],
                     DatamapValue[1]},
               },
               {
                  DataMinValue=><-.DataRange[0].DataMaxValue,,,,,DataMaxValue=0.,,,,,,,,,,,,controlPoints=>
                  {DatamapValue[2],
                     DatamapValue[3]},
               },
               {
                  DataMinValue=><-.DataRange[1].DataMaxValue,,,,,DataMaxValue=1.,,,,,,,,,,,,controlPoints=>
                  {DatamapValue[4],
                     DatamapValue[5]},
               }};
            DatamapValue[6] = {
               
               {
                  v1=0.,v2=0.6348015666,v3=0.8271973729,v4=1.,value=-1.
               },
               {
                  v1=0.,v2=0.6666666865,v3=1.,v4=1.,value=0.
               },
               {
                  v1=0.,v2=0.6666666865,v3=1.,v4=1.,value=0.
               },
               {
                  v1=1.,v2=0.,v3=1.,v4=1.,value=0.
               },
               {
                  v1=1.,v2=0.,v3=1.,v4=1.,value=0.
               },
               {
                  v1=0.9887584448,v2=0.07732812315,v3=0.3479741216,v4=1.,value=0.
               }
            };
         };
         regen = 0;
         delta = 0.007843137719;
      };
      Presets {
         preset_load {
            preset = {
               0.,0.,0.,0.,1.,0.5,0.,0.,0.,1.,0.5,1.,1.,0.,0.,1.,1.,1.,
0.,0.
            };
         };
      };
   };
   CMAP_EDTR.ColorMapEditor fluid_velocity<NEx=72.,NEy=459.> {
      InObj => <-.glyph.out_obj;
      TransferFunction {
         PointEditor {
            real_scalar = 2.260397923e-05;
         };
         Scene {
            Top {
               AltObject<instanced=0>;
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      UIshell {
         x = 1110;
         y = 110;
         iconic = 1;
      };
      DmapRamp {
         Scene {
            Top {
               AltObject<instanced=0>;
               Xform {
                  ocenter = {170.,16.,0.};
                  dcenter = {-3122.,-267.13,0.};
               };
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      EditFields {
         UIoptionBoxRefresh {
            selectedItems = {0};
         };
         UIoptionRefresh {
            set = 1;
         };
      };
      ColorWheel {
         Scene {
            Top {
               AltObject<instanced=0>;
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
         HSPosConverter {
            hue = 0.7097604324;
            saturation = 0.04079174995;
         };
         PickField {
            xform {
               xlate = {-0.010204,-0.0394949,
0.};
            };
         };
      };
      ValueCone {
         Scene {
            Top {
               AltObject<instanced=0>;
               Xform {
                  ocenter = {0.15,0.51,0.};
                  dcenter = {26.5619,90.31,
0.};
               };
            };
            View {
               View {
                  trigger = 1;
               };
            };
         };
      };
      MakeDmap {
         dmap_out {
            dataMax = 0.;
            DataRange[2] = {
               {
                  DataMinValue=2.260397923e-05,,,,,DataMaxValue=0.002259901259,,,,,,,,,,,,controlPoints=>
                  {DatamapValue[0],
                     DatamapValue[1]},
               },
               {
                  DataMinValue=><-.DataRange[0].DataMaxValue,,,,,DataMaxValue=0.02239557914,,,,,,,,,,,,controlPoints=>
                  {DatamapValue[2],
                     DatamapValue[3]},
               }};
            DatamapValue[4] = {
               {
                  v1=1.,v2=0.7097603083,v3=0.04079173505,v4=1.,value=2.260397923e-05
               },
               {
                  v1=0.,v2=0.,v3=1.,v4=1.,value=0.002259901259
               },
               {
                  v1=0.,v2=0.,v3=1.,v4=1.,value=0.002259901259
               },
               {
                  v1=0.,v2=0.,v3=1.,v4=1.,value=0.02239557914
               }};
         };
         regen = 0;
      };
      Presets {
         preset_load {
            preset = {
               0.,1.,0.,0.,1.,0.1,0.,1.,0.,0.,1.,0.,1.,0.,0.
            };
         };
      };
   };
};
