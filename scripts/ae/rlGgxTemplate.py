import pymel.core as pm
import mtoa.utils as utils
import mtoa.ui.ae.utils as aeUtils
from mtoa.ui.ae.shaderTemplate import ShaderAETemplate

class AErlGgxTemplate(ShaderAETemplate):

    def setup(self):
        # Add the shader swatch to the AE
        self.addSwatch()
        self.beginScrollLayout()

        # Add a list that allows to replace the shader for other one
        # self.addCustom('message', 'AEshaderTypeNew', 'AEshaderTypeReplace')

        # Begins a "Color Section"
        self.beginLayout("Common Attributes", collapse=False)
        self.addControl("roughness");
        self.addControl("ior", label="IOR");
        self.addControl("opacity");
        self.endLayout()

        self.beginLayout("Diffuse", collapse=False)
        self.addControl("KdColor", label="Color");
        self.addControl("Kd", label="Weight");
        self.endLayout()

        self.beginLayout("Specular", collapse=False)
        self.addControl("KsColor", label="Color");
        self.addControl("Ks", label="Weight");
        self.endLayout()

        self.beginLayout("Refract", collapse=False)
        self.addControl("KtColor", label="Color");
        self.addControl("Kt", label="Weight");
        self.endLayout()

        # include/call base class/node attributes
        pm.mel.AEdependNodeTemplate(self.nodeName)

        # Add Section for the extra controls not displayed before
        self.addExtraControls()
        self.endScrollLayout()
