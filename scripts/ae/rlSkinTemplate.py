import pymel.core as pm
import mtoa.utils as utils
import mtoa.ui.ae.utils as aeUtils
from mtoa.ui.ae.shaderTemplate import ShaderAETemplate

class AErlSkinTemplate(ShaderAETemplate):

    def setup(self):
        # Add the shader swatch to the AE
        self.addSwatch()
        self.beginScrollLayout()

        # Add a list that allows to replace the shader for other one
        # self.addCustom("message", "AEshaderTypeNew", "AEshaderTypeReplace")

        # Begins a "Color Section"
        self.beginLayout("Sheen", collapse=False)
        self.addControl("sheen_color", label="Color")
        self.addControl("sheen_weight", label="Weight")
        self.addControl("sheen_roughness", label="Roughness")
        self.addControl("sheen_ior", label="IOR")
        self.endLayout()

        self.beginLayout("Specular", collapse=False)
        self.addControl("specular_color", label="Color")
        self.addControl("specular_weight", label="Weight")
        self.addControl("specular_roughness", label="Roughness")
        self.addControl("specular_ior", label="IOR")
        self.endLayout()

        self.beginLayout("SSS", collapse=False)
        self.addControl("sss_color", label="Color")
        self.addControl("sss_weight", label="Weight")
        self.addControl("sss_dist_multiplier", label="Distance Multiplier")
        self.addControl("sss_scatter_dist", label="Scatter Distance")
        self.addControl("sss_cavity_fadeout", label="Cavity Affects Diffusion")
        self.endLayout()

        self.beginLayout("Opacity", collapse=True)
        self.addControl("opacity", label="Weight")
        self.addControl("opacity_color", label="Color")
        self.endLayout()

        self.addBumpLayout()
        self.addAOVLayout(aovReorder=["sheen", "specular", "sss"])

        # include/call base class/node attributes
        pm.mel.AEdependNodeTemplate(self.nodeName)

        # Add Section for the extra controls not displayed before
        self.addExtraControls()
        self.endScrollLayout()
