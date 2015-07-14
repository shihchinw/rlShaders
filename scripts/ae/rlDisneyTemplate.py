import pymel.core as pm
import mtoa.utils as utils
import mtoa.ui.ae.utils as aeUtils
from mtoa.ui.ae.shaderTemplate import ShaderAETemplate

class AErlDisneyTemplate(ShaderAETemplate):

    def setup(self):
        # Add the shader swatch to the AE
        self.addSwatch()
        self.beginScrollLayout()

        # Add a list that allows to replace the shader for other one
        # self.addCustom('message', 'AEshaderTypeNew', 'AEshaderTypeReplace')

        # Begins a "Color Section"
        self.beginLayout("Material Attributes", collapse=False)
        self.addControl("base_color");
        self.addControl("subsurface");
        self.addControl("metallic");
        self.addControl("specular");
        self.addControl("specular_tint");
        self.addControl("roughness");
        self.addControl("anisotropic");
        self.addControl("sheen");
        self.addControl("sheen_tint");
        self.addControl("clearcoat");
        self.addControl("clearcoat_gloss");
        self.endLayout()

        self.addBumpLayout()

        # include/call base class/node attributes
        pm.mel.AEdependNodeTemplate(self.nodeName)

        # Add Section for the extra controls not displayed before
        self.addExtraControls()
        self.endScrollLayout()
