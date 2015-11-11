import pymel.core as pm
from mtoa.ui.ae.shaderTemplate import ShaderAETemplate


class AErlDisneyTemplate(ShaderAETemplate):

    def setup(self):
        # Add the shader swatch to the AE
        self.addSwatch()
        self.beginScrollLayout()

        self.beginLayout("Material Attributes", collapse=False)
        self.addControl("base_color")
        self.addControl("subsurface")
        self.addControl("metallic")
        self.addControl("specular")
        self.addControl("specular_tint")
        self.addControl("roughness")
        self.addControl("anisotropic")
        self.addControl("sheen")
        self.addControl("sheen_tint")
        self.addControl("clearcoat")
        self.addControl("clearcoat_gloss")
        self.addControl("opacity")
        self.endLayout()

        self.beginLayout("Extended Controls", collapse=False)
        self.addControl("indirectDiffuseScale")
        self.addControl("indirectSpecularScale")
        self.endLayout()

        self.addBumpLayout()
        self.addAOVLayout(aovReorder=["direct_diffuse", "direct_specular",
                          "indirect_diffuse", "indirect_specular"])

        # include/call base class/node attributes
        pm.mel.AEdependNodeTemplate(self.nodeName)

        # Add Section for the extra controls not displayed before
        self.addExtraControls()
        self.endScrollLayout()
